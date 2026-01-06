classdef radarClass < handle
    properties
        %instance details
        ID
        sceneario
        param_Q=250; %OLD WAS 55
        param_R=150;       % OLD WAS 150.5
        param_P=100; % OLD WAS 5
        %basic signals
        radar_i
        radar_q
        radar_dist
        ecg_gt
        resp_gt
        %statistics
        HrEstVar
        HrEstRoughNorm
        IQRadCV
        DistSNR_HR_dB
        DistPeakToFloor
        DistHrPeakHz

        %proccessed signals
        radar_decimated
        HrSignal
        KF_HrSignal
        RrSignal
        HrPeaks
        HrPeaksAfterKalman
        RrPeaks
        ecgPeaks
        Rrpeaks_gt

        %results
        HrEst
        HrEstAfterMedian
        HrGtEstAfterMedian
        HrEstAfterKalman
        mechanicalDelayKalman
        mechanicalDelayMedian 
        %correlating by timestamp
        CorrKalmanHr_on_gt_time %kalman at gt timeline
        CorrMedianHr_on_gt_time %Median HR at gt timeline
        CorrGt %gt at gt timeline
        kalmanCorrValue
        medianCorrValue
        %TODO- show results with these 2 vectors
        RrEst
        HrGtEst
        RrGtEst
        HrGtMean
        HrGtDiff
        mseRaw
        mseFitted
        maeRaw
        maeFitted
        mse2HrRaw
        mse2Hr
        mae2Hr
        rawErr2Hr
        corrtime


        %additional information
        vTimeOriginal
        vTimeNew
        
        fs_radar
        fs_new
        fs_ecg
        % filters?
        
    end
    properties(Constant)
        %MaxDiffFromGt = 0.4; %maximum difference between beats of radar and GT allowed 

    end
    methods
        % constructor: call when making a new object
        function obj = radarClass(ID,scenario,fs,ecg,radar_i, radar_q,gt_resp)
            arguments
                ID            
                scenario (1,1) string {mustBeMember(scenario, ["Resting","Valsalva","Apnea","TiltDown","TiltUp"])} 
                fs {mustBeNonnegative}
                ecg {mustBeColumn}
                radar_i {mustBeColumn} 
                radar_q {mustBeColumn} = 0 %optional, in case we get radar_dist
                gt_resp{mustBeColumn} = 0 
            end
                 %TODO : create the process for getting TFM_respiration

            %similliar to getvitalsigns() instead of getting it directly
            obj.ID = string(ID);
            obj.fs_radar = fs;
            obj.fs_ecg = fs;
            obj.sceneario = scenario;
            obj.ecg_gt = ecg;
            obj.resp_gt = gt_resp;
            if(radar_q==0)
                obj.radar_dist=radar_i;
                sprintf('Only one signal detected, saved as radar_dist')
            else
            obj.radar_i = radar_i;
            obj.radar_q = radar_q;
            end
            obj.vTimeOriginal= 1/fs:1/fs:(length(ecg))/fs; % len-1?

        end
        
        %I-Q compensation- create radar_dist from radar_i and radar_Q
            %TODO

        % downSample radar dist- accept new fs or use default. change
        % fs_radar accordingly
        function DS = DownSampleRadar(obj,fs)
            if(nargin<2 || isempty(fs))
                fs=100;

            end
            DS=obj.fs_radar/fs;
            %sprintf('Set new radar fs to %d',fs)
            obj.radar_decimated = decimate(obj.radar_dist, DS);
            obj.fs_new = fs;
            obj.vTimeNew = 1/fs:1/fs:length(obj.radar_decimated)/fs; %len-1?
        end
    %end 
    %methods(Static)
        % apply HR filters and make %TODO use outside filter so we wont
        % designfilt each iteration
        function HrFilter(obj, LPF, HPF)
            if nargin < 3
                 HPF = [];
            end
            if nargin < 2
                LPF = [];
            end
            if((obj.fs_new)==0)
                obj.fs_new=obj.fs_radar;
                sprintf('warning: signal not yet decimated')
            end
            if(isempty(HPF) || nargin<3)
                firH = designfilt('highpassfir','StopbandFrequency',0.4,...
                'PassbandFrequency',0.7,'StopbandAttenuation',60, ...
                'SampleRate',obj.fs_new);
            else
                firH=HPF;
            end
            if(isempty(LPF) || nargin<2)
                firL = designfilt('lowpassfir','PassbandFrequency',3,...
                'StopbandFrequency',4,'StopbandAttenuation',80, ...
                'SampleRate',obj.fs_new);
            else
                firL=LPF;
            end
            % Apply the filters to the decimated radar signal
            obj.HrSignal = filtfilt(firL, obj.radar_decimated);
            obj.HrSignal = filtfilt(firH, obj.HrSignal);
        end

        %
        function RrFilter(obj, LPF, HPF) %currently without highpass. maybe pass through median filter. %TODO use outside filter so we wont
        % designfilt each iteration
          if nargin < 2
                LPF = [];
          end
          if(isnan(obj.fs_new))
              obj.fs_new=obj.fs_radar;
              sprintf('warning: signal not yet decimated')
          end
          
            % firH = designfilt('highpassfir','StopbandFrequency',0.02,...
            % 'PassbandFrequency',0.1,'StopbandAttenuation',40, ...
            % 'SampleRate',obj.fs_new);
            if(isempty(LPF) || nargin<2)
                firL = designfilt('lowpassfir','PassbandFrequency',0.4,...
                'StopbandFrequency',0.8,'StopbandAttenuation',40, ...
                'SampleRate',obj.fs_new);
            else
            firL=LPF;
            firH=HPF;
            end
            % Apply the filters to the decimated radar signal
            obj.RrSignal = filtfilt(firL, obj.radar_decimated);
            obj.RrSignal = filtfilt(firH, obj.RrSignal);
        end
        %finding the peaks and normalize to seconds
        function FindPeaks(obj)
            thresholdHr= mean(abs((obj.HrSignal)))*0.05;
            thresholdRr= mean(abs((obj.RrSignal)))*0.05;
            [~,obj.HrPeaks, ~,~] = findpeaks(obj.HrSignal, "MinPeakHeight",...
        thresholdHr,'MinPeakDistance',0.33*obj.fs_new);
            [~,obj.RrPeaks, ~,~] = findpeaks(obj.RrSignal, "MinPeakHeight",...
        thresholdRr,'MinPeakDistance',2*obj.fs_new);
            [~,obj.ecgPeaks,~] = pan_tompkin(obj.ecg_gt,obj.fs_ecg,0); 
            [~,obj.Rrpeaks_gt, ~,~] = findpeaks(obj.resp_gt, "MinPeakHeight",...
        thresholdRr,'MinPeakDistance',2*(obj.fs_new/2.5));          
            obj.HrPeaks = obj.HrPeaks / obj.fs_new;
            obj.RrPeaks = obj.RrPeaks / obj.fs_new; %in seconds 
            obj.ecgPeaks = obj.ecgPeaks / obj.fs_ecg; % in seconds 
            obj.Rrpeaks_gt = obj.Rrpeaks_gt /(100);
        end

        % finding the rates
        function FindRates(obj)
            obj.HrEst = 60 ./  diff(obj.HrPeaks);
            obj.HrGtEst = 60 ./  diff(obj.ecgPeaks); 
            obj.RrEst = 60 ./  diff(obj.RrPeaks);
            obj.RrGtEst = 60 ./ diff(obj.Rrpeaks_gt);            
        end
        function ComputePreFilterStats(obj)
            % ComputePreFilterStats
            % Computes statistics available BEFORE running the Kalman filter.
            % Stores results in:
            %   HrEstVar, HrEstRoughNorm, IQRadCV, DistHrPeakHz, DistSNR_HR_dB, DistPeakToFloor
    
            % ---------- HR estimate stats (from HrEst) ----------
            obj.HrEstVar = NaN;
            obj.HrEstRoughNorm = NaN;
            
            if ~isempty(obj.HrEst)
                hr = obj.HrEst(:);
                hr = hr(~isnan(hr) & isfinite(hr));
                if numel(hr) >= 5
                    obj.HrEstVar = var(hr, 1); % population variance (stable for comparisons)
                    rough = mean(abs(diff(hr)));
                    obj.HrEstRoughNorm = rough / max(mean(hr), eps);
                end
            end
            
            % ---------- Radar IQ "radial" dispersion stats ----------
            % Uses raw radar_i/q if present, otherwise radar_decimated if it's complex
            obj.IQRadCV  = NaN;
            
            
            [ri, rq] = getIQ(obj);
            if ~isempty(ri) && ~isempty(rq)
                ri = ri(:); rq = rq(:);
                n = min(numel(ri), numel(rq));
                ri = ri(1:n); rq = rq(1:n);
            
                % Robust center
                ci = median(ri, 'omitnan');
                cq = median(rq, 'omitnan');
            
                rad = hypot(ri - ci, rq - cq);
                rad = rad(~isnan(rad) & isfinite(rad));
            
                if numel(rad) >= 20
                    medr = median(rad);
                    iqr_r = iqr(rad);
                    obj.IQRadCV = iqr_r / max(medr, eps);
            
                    % If you also want MAD:
                    % obj.IQRadMAD = mad(rad, 1); % 1 => normalized by median
                end
            end
            
            % ---------- Spectral "dist" metrics on HrSignal ----------
            % DistHrPeakHz: peak frequency (Hz) in HR band
            % DistSNR_HR_dB: simple peak-vs-floor dB measure
            % DistPeakToFloor: linear peak/floor
            obj.DistHrPeakHz   = NaN;
            obj.DistSNR_HR_dB  = NaN;
            obj.DistPeakToFloor = NaN;
            
            if ~isempty(obj.HrSignal) && ~isempty(obj.fs_new) && obj.fs_new > 0
                x = obj.HrSignal(:);
                x = x(~isnan(x) & isfinite(x));
            
                if numel(x) >= 128
                    fs = obj.fs_new;
            
                    % Detrend and window
                    x = x - mean(x);
                    w = hann(numel(x));
                    X = abs(fft(x .* w));
                    f = (0:numel(X)-1).' * (fs/numel(X));
            
                    % HR band (Hz): 0.7–3.0 ~ 42–180 bpm (adjust if you want)
                    f1 = 0.7; f2 = 3.0;
                    band = (f >= f1) & (f <= f2);
            
                    if any(band)
                        Xb = X(band);
                        fb = f(band);
            
                        [pk, idxPk] = max(Xb);
                        obj.DistHrPeakHz = fb(idxPk);
            
                        % Floor estimate: median in band (robust)
                        floorVal = median(Xb, 'omitnan');
            
                        obj.DistPeakToFloor = pk / max(floorVal, eps);
                        obj.DistSNR_HR_dB   = 20*log10(obj.DistPeakToFloor); % amplitude ratio -> dB
                    end
                end
            end
    end

        % ---- helper inside class (put under methods(Access=private) if you prefer) ----
        function [ri, rq] = getIQ(obj)
            ri = [];
            rq = [];
            if ~isempty(obj.radar_i) && ~isempty(obj.radar_q)
                ri = obj.radar_i;
                rq = obj.radar_q;
                return;
            end
            % If radar_decimated is complex, split it
            if ~isempty(obj.radar_decimated) && ~isreal(obj.radar_decimated)
                z = obj.radar_decimated;
                ri = real(z);
                rq = imag(z);
            end
        end

        function [Q,R] = ProduceKalmanCoeff(obj)
            variance = var(obj.HrEst);
            
            variance = max(variance, 1e-4);   % safety floor
            scenario = obj.sceneario;
            R= 19.31;
            switch string(scenario)
                case "Resting"
                    C = 1.0e3;  a = -1.0;
                    R= 9.21;
                case "Apnea"
                    C = 1.5e3;  a = -1.2;
                    R= 3.43;

                case "TiltDown"
                    C = 1.0e3;  a = -0.7;
                    R= 11.79;

                case "TiltUp"
                    C = 2.0e3;  a = -1.4;
                    R= 35;

                case "Valsalva"
                    C = 2.5e3;  a = -1.5;
                    R = R * variance.^(-0.3);


                otherwise
                    C = 1.2e3;  a = -1.0;  % fallback
            end

            Q = C * variance.^a;


        % Optional bounds for numerical stability
            Q = min(max(Q, 1e-3), 3e3);
        end



        function MedianHr(obj,size)
            if(nargin<2)
                size=10;
            end
          gtHr=obj.HrGtEst(:);
          calculatedHr=obj.HrEst(:);
          %minlen=min(length(gtHr),length(calculatedHr));

          gt_after_median_filt=medfilt1(gtHr,size);
          calculatedHr_after_median_filt=medfilt1(calculatedHr,size);

          %FOR RADAR SIGNAL
          upper_bound=1.4; %was 1.2
          lower_bound=0.6; %was 0.8
         indx_radar = find(calculatedHr > upper_bound * calculatedHr_after_median_filt | calculatedHr < lower_bound * calculatedHr_after_median_filt); 
         indx_gt = find(gtHr > upper_bound * gt_after_median_filt | gtHr < lower_bound * gt_after_median_filt); 
                    
         hr_replaced_with_median_in_outliers = obj.replaceOutliersByInterp(calculatedHr, indx_radar);
         gt_replaced_with_median_in_outliers = obj.replaceOutliersByInterp(gtHr, indx_gt);
         hr_replaced_with_median_in_outliers=hr_replaced_with_median_in_outliers(:);
         gt_replaced_with_median_in_outliers=gt_replaced_with_median_in_outliers(:);
        % cor_after_median_Filter_on_both_signals = corr(hr_replaced_with_median_in_outliers,gt_replaced_with_median_in_outliers)
        obj.HrEstAfterMedian = hr_replaced_with_median_in_outliers;
        obj.HrGtEstAfterMedian = gt_replaced_with_median_in_outliers;
        end
        function [] = CalcError(obj,hrToCompare)
            min_len = min(length(obj.CorrGt), length(hrToCompare));
            v1=hrToCompare(1:min_len);
            v2=obj.HrGtEstAfterMedian(1:min_len);


            v2=v2(:);
            mseraw= (v1-v2).^2;
            obj.mseRaw = mean(mseraw);
            obj.maeRaw = mean(abs(v1-v2));
            obj.mse2HrRaw = sortrows([v2, mseraw]); % N,2, acending error per beat rate
            gt_vs_mse_raw_matrix = [v2, mseraw];
            gt_vs_rawerr_raw_matrix = [v2, v1-v2];
            gt_vs_abserr_raw_matrix = [v2, abs(v1-v2)];
            obj.mse2Hr = sortrows(gt_vs_mse_raw_matrix); % N,2, acending error per beat rate
            obj.rawErr2Hr = sortrows(gt_vs_rawerr_raw_matrix ); 
            obj.mae2Hr = sortrows(gt_vs_abserr_raw_matrix );   
        end


function timeFitting(obj)
    tEstPk = obj.HrPeaks(:);                 % seconds
    hrEst  = obj.HrEstAfterKalman;             % bpm
    tEstHr = (tEstPk(1:end-1) + tEstPk(2:end))/2;  % midpoint time per IBI
    tGtPk = obj.ecgPeaks(:);                % seconds
    hrGt  = obj.HrGtEstAfterMedian;
    %hrGt  = obj.HrGtEst(:);
    tGtHr = (tGtPk(1:end-1) + tGtPk(2:end))/2;
        tMedian = obj.HrEstAfterMedian;
      % Overlapping window
    t0 = max(tEstHr(1), tGtHr(1));
    t1 = min(tEstHr(end), tGtHr(end));
    
    idxGt = (tGtHr >= t0) & (tGtHr <= t1);
    
    tGrid = tGtHr(idxGt);   % GT timeline
    obj.corrtime = tGrid;
    hrEst_on_GT = interp1(tEstHr, hrEst, tGrid, 'pchip');
    hrMed_on_GT = interp1(tEstHr, tMedian, tGrid, 'pchip');
    hrGt_on_GT  = hrGt(idxGt);
    obj.CorrKalmanHr_on_gt_time=hrEst_on_GT;
    obj.CorrGt= hrGt_on_GT;
    obj.CorrMedianHr_on_gt_time = hrMed_on_GT;
    obj.kalmanCorrValue = corr(hrEst_on_GT, hrGt_on_GT, 'Rows','complete');
    obj.medianCorrValue = corr(hrMed_on_GT, hrGt_on_GT, 'Rows','complete');

    fprintf('  CORRELATION ANALYSIS: %s - %s\n',obj.ID, obj.sceneario);
    fprintf('  %-25s : %.4f\n', 'Kalman Filter Correlation', obj.kalmanCorrValue);
    fprintf('  %-25s : %.4f\n', 'Median Filter Correlation', obj.medianCorrValue);
    fprintf('  %-25s : %d\n', 'Total Valid Samples', sum(isfinite(hrEst_on_GT)));   
end

%% KALMAN 
function KalmanFilterBeats(obj,Q,R_base)
    % Improved Kalman Filter for Heart Rate
    % 1. Robust against NaN/Empty inputs.
    % 2. Uses Phase-Locking to prevent drift.
    % 3. Filters Heart Rate (BPM) while preserving beat alignment.

    % ---- 1. Guards: Check for valid input ----
    % We need at least 2 peaks to calculate an interval
    if isempty(obj.HrPeaks) || length(obj.HrPeaks) < 2
        obj.HrPeaksAfterKalman = obj.HrPeaks;
        obj.HrEstAfterKalman = []; 
        warning('Kalman:InputSparse', 'Input HrPeaks too sparse (<2). Returning original.');
        return;
    end

    fs = obj.fs_new;
    if isempty(fs) || fs == 0, fs = 100; end % Safety fallback

    % Prepare Input Data
    raw_pk_times = unique(sort(obj.HrPeaks(:)), 'stable');
    ibi_meas = diff(raw_pk_times);
    hr_meas = 60 ./ ibi_meas; 

    % ---- 2. Kalman Initialization (State: [HR]) ----
    % We use a 1-state Random Walk model. Velocity state (HRdot) 
    % often overshoots on radar data due to sudden noise.

    % Find a valid starting HR (median of first few valid beats)
    valid_mask = (hr_meas > 40 & hr_meas < 200 & isfinite(hr_meas));
    valid_data = hr_meas(valid_mask);

    if isempty(valid_data)
        x = 75; % Fallback if signal is 100% garbage
    else
        x = median(valid_data(1:min(5, end)));
    end
    P = 65.0; 
    if(nargin<3)
                  % Initial uncertainty 10
        Q = 35;          % Process noise (Standard deviation of beat-to-beat change) oldval 0.5
        R_base = 112;     % Measurement noise (Trust in the radar peak location) oldval 5 
        
    end
    
    hr_hat = nan(size(hr_meas));

    % ---- 3. Filtering Loop ----
    for k = 1:length(hr_meas)
        z = hr_meas(k);

        % Predict
        x_pred = x;
        P_pred = P + Q;

        % Validate Measurement
        % 1. Is it physically possible? (40-200 BPM)
        % 2. Is it within the statistical gate? (Mahalanobis distance)
        is_physio = (z >= 40 && z <= 200 && isfinite(z));

        innovation = z - x_pred;
        S = P_pred + R_base;

        % Gate: If error > 3 standard deviations, it's an outlier (or missed beat)
        if (innovation^2) > (9 * S) 
            gate_passed = false;
        else
            gate_passed = true;
        end
        % Update Step
        if is_physio && gate_passed
            K = P_pred / S;
            x = x_pred + K * innovation;
            P = (1 - K) * P_pred;
        else
            % Reject measurement, trust model history
            x = x_pred; 
            P = P_pred; 
        end

        hr_hat(k) = x;
    end

    obj.HrEstAfterKalman = hr_hat;

    % ---- 4. Phase-Locked Reconstruction ----
    % Prevents the "drifting away" issue of simple cumsum.
    % We predict where the next beat *should* be. If a raw peak is close,
    % we snap to it. If not, we use the prediction (smoothing).

    rec_peaks = zeros(length(hr_hat)+1, 1);
    rec_peaks(1) = raw_pk_times(1); % Anchor to first peak
    curr_time = raw_pk_times(1);

    for k = 1:length(hr_hat)
        % Predict IBI based on filtered HR
        pred_ibi = 60 / hr_hat(k);
        next_pred = curr_time + pred_ibi;

        % Check if a raw peak exists near this prediction
        % Window: +/- 20% of the beat interval
        window = 0.2 * pred_ibi;
        [min_err, idx] = min(abs(raw_pk_times - next_pred));

        if min_err < window
            % Found a matching raw peak -> Snap to it (Corrects Phase)
            next_time = raw_pk_times(idx);
        else
            % No match -> Use Kalman prediction (Smooths over noise/misses)
            next_time = next_pred;
        end
        rec_peaks(k+1) = next_time;
        curr_time = next_time;
    end

    obj.HrPeaksAfterKalman = rec_peaks;
end
function out = kalmanSmoothRadarDist(obj, config)
% kalmanSmoothRadarDist
% Kalman smoothing for raw radar_dist at high fs (e.g. 2000 Hz).
% Supports 2-state (x, xdot) and 3-state (x, xdot, xddot).
%
% INPUTS
%   radar_dist : Nx1 (or 1xN) double, raw distance (meters or arbitrary units)
%   fs         : sampling rate in Hz (e.g. 2000)
%   config     : struct with optional fields:
%       .order            (2 or 3) default 3
%       .useRTS           (true/false) default true (Rauch-Tung-Striebel smoother)
%       .initPScale       default 10
%       .sigmaMeas        measurement noise std (same units as radar_dist)
%       .sigmaAcc         process accel noise std (units/s^2) for order=3
%       .sigmaJerk        process jerk noise std (units/s^3) for order=2
%       .robustDetrend    (true/false) default true
%       .hpCutHz          detrend highpass cutoff (Hz) default 0.05
%       .nanPolicy        'interp' or 'omit' default 'interp'
%
% OUTPUT (struct)
%   out.x, out.xdot, out.xddot (if order=3)
%   out.x_filt (filtered forward pass), out.x_smooth (RTS)
%   out.innov, out.S (innovation & its variance)
%   out.params (A,H,Q,R,dt,order, estimated sigmas)

arguments
    obj {mustBeNonempty}
    config.order (1,1) double {mustBeMember(config.order,[2 3])} = 3
    config.useRTS (1,1) logical = true
    config.initPScale (1,1) double {mustBePositive} = 10
    config.sigmaMeas (1,1) double = NaN
    config.sigmaAcc (1,1) double = NaN
    config.sigmaJerk (1,1) double = NaN
    config.robustDetrend (1,1) logical = true
    config.hpCutHz (1,1) double {mustBeNonnegative} = 0.05
    config.nanPolicy (1,:) char {mustBeMember(config.nanPolicy,{'interp','omit'})} = 'interp'
end


x = obj.HrSignal(:);
fs = obj.fs_new;
N = numel(x);
dt = 1/fs;

% --- NaN handling ---
nanMask = isnan(x);
if any(nanMask)
    switch config.nanPolicy
        case 'interp'
            t = (0:N-1)'*dt;
            x(nanMask) = interp1(t(~nanMask), x(~nanMask), t(nanMask), 'linear', 'extrap');
        case 'omit'
            % For 'omit', we keep NaNs and skip measurement update when NaN
            % (implemented below)
    end
end

% --- Robust detrend / drift removal (VERY important at 2000 Hz) ---
x0 = x;
if config.robustDetrend
    % remove slow drift with a gentle highpass implemented via moving median
    % window length ~ 1/config.hpCutHz seconds (clamped)
    wSec = max(5, 1/max(config.hpCutHz,1e-3)); % seconds
    w = max(201, 2*floor((wSec*fs)/2)+1);      % odd, >=201
    trend = movmedian(x, w, 'omitnan');
    x = x - trend;
end

% --- Estimate measurement noise if not provided ---
% Use robust estimate on high-frequency component: diff(x) is dominated by noise at high fs.
if isnan(config.sigmaMeas)
    % Estimate measurement noise from short-window residual (robust)
    L = max(5, round(0.05*fs));          % ~50 ms window (fs dependent)
    x_lp = movmean(x, L, 'omitnan');
    resid = x - x_lp;
    sigmaMeas = mad(resid,1);
    if sigmaMeas == 0 || ~isfinite(sigmaMeas)
        sigmaMeas = std(resid) + eps;
    end
else
    sigmaMeas = config.sigmaMeas;
end
sigmaMeas = max(sigmaMeas, 1e-4*std(x) + eps);
sigScale = std(x(~isnan(x)));
sigmaFloor = 0.5 * sigScale;      % 1% of signal std (start here)
sigmaMeas  = max(sigmaMeas, sigmaFloor);
R = sigmaMeas^2;




% --- Build state-space model ---
order = config.order;
H = zeros(1, order); H(1) = 1; % measure position only

if order == 2
    % State: [x; xdot]
    A = [1 dt;
         0 1];

    % Process noise driven by "jerk" (random accel changes) on xdot
    % Continuous-time white noise on acceleration integrated gives Q ~ sigmaJerk^2 * [...]
    if isnan(config.sigmaJerk)
        % crude default: set jerk so that model is moderately flexible
        sigmaJerk = 10*sigmaMeas / (dt^2); % heuristic
    else
        sigmaJerk = config.sigmaJerk;
    end
    q = sigmaJerk^2;
    Q = q * [dt^3/3, dt^2/2;
             dt^2/2, dt];

elseif order == 3
    % State: [x; xdot; xddot]
    A = [1 dt 0.5*dt^2;
         0 1  dt;
         0 0  1];

    % Process noise driven by "acceleration random walk" (white noise on jerk)
    if isnan(config.sigmaAcc)
        % default: allow moderate accel changes; scale from sigmaMeas
        sigmaAcc = 5*sigmaMeas / (dt^2); % heuristic (units/s^2)
    else
        sigmaAcc = config.sigmaAcc;
    end
    q = sigmaAcc^2;
    % Standard discrete Q for constant-acceleration model with white noise on acceleration
    Q = q * [dt^5/20, dt^4/8,  dt^3/6;
             dt^4/8,  dt^3/3,  dt^2/2;
             dt^3/6,  dt^2/2,  dt];
end

R = sigmaMeas^2;

% --- Initialize ---
xhat = zeros(order, N);
P = eye(order) * config.initPScale;

% Initialize state from first samples
xhat(1,1) = x(1);
if N >= 2
    xhat(2,1) = (x(2)-x(1))/dt; % crude derivative init
end
if order == 3 && N >= 3
    xhat(3,1) = ((x(3)-x(2)) - (x(2)-x(1))) / (dt^2);
end

innov = nan(1,N);
S = nan(1,N);

% Store forward pass for RTS
xpred_store = zeros(order,N);
Ppred_store = zeros(order,order,N);
P_store     = zeros(order,order,N);

% --- Forward Kalman filter ---
for k = 2:N
    % Predict
    xpred = A * xhat(:,k-1);
    Ppred = A * P * A' + Q;

    xpred_store(:,k) = xpred;
    Ppred_store(:,:,k) = Ppred;

    zk = x0(k);
    hasMeas = ~(isnan(zk) && strcmp(config.nanPolicy,'omit'));

    if hasMeas
        % Use detrended measurement if detrending enabled, else raw
        z = x(k);

        % Update
        y = z - H*xpred;             % innovation
        Sk = H*Ppred*H' + R;         % innovation variance (scalar)
        K = (Ppred*H') / Sk;         % Kalman gain

        xnew = xpred + K*y;
        Pnew = (eye(order) - K*H) * Ppred;

        innov(k) = y;
        S(k) = Sk;
    else
        xnew = xpred;
        Pnew = Ppred;
    end

    xhat(:,k) = xnew;
    P = Pnew;
    P_store(:,:,k) = P;
    if k == 2
        fprintf('sigmaMeas=%.3g, R=%.3g\n', sigmaMeas, R);
        fprintf('Sk=%.3g, Kpos=%.3g\n', Sk, K(1));
    end
end

% --- RTS smoother (optional) ---
xs = xhat;
Ps = P_store;
if config.useRTS
    for k = N-1:-1:1
        Pk = P_store(:,:,k);
        Ppred_next = Ppred_store(:,:,k+1);
        if all(Ppred_next(:)==0)
            continue
        end
        G = (Pk*A') / Ppred_next; % smoother gain
        xs(:,k) = xhat(:,k) + G*(xs(:,k+1) - xpred_store(:,k+1));
        Ps(:,:,k) = Pk + G*(Ps(:,:,k+1) - Ppred_next)*G';
    end
end

% --- Outputs ---
out = struct();
out.params = struct('A',A,'H',H,'Q',Q,'R',R,'dt',dt,'fs',fs,'order',order, ...
                    'sigmaMeas',sigmaMeas);

out.x_filt = xhat(1,:).';
out.x_smooth = xs(1,:).';
out.innov = innov(:);
out.S = S(:);

if order >= 2
    out.xdot_filt = xhat(2,:).';
    out.xdot_smooth = xs(2,:).';
end
if order == 3
    out.xddot_filt = xhat(3,:).';
    out.xddot_smooth = xs(3,:).';
end
obj.KF_HrSignal = out.x_smooth(:);
end

function [delay_sec_normal,sign]= FindMechDelay(obj)
    %first, the simplest method - XCORR on the HR peaks times of different
    %signals and ECG
    hr_peak_times  = obj.HrPeaks; %sample rate - 100 
    %hr_peak_times_after_kalman = obj.HrPeaksAfterKalman; %sample rate - 100 

    gt_peak_times = obj.ecgPeaks ; % sample rate  - 2000

    fs_common = 1000;
    max_time = max([max(hr_peak_times), max(gt_peak_times)]);
   % max_time_betweenKalman = max([max(hr_peak_times_after_kalman), max(gt_peak_times)]);


    grid_len_peaks_gt = ceil(max_time*fs_common)+1;
    %grid_len_peaksAfterKalman_gt =  ceil(max_time_betweenKalman*fs_common)+1;

    sig_ECG = zeros(grid_len_peaks_gt, 1);
    sig_HR  = zeros(grid_len_peaks_gt, 1);

    %sig_ECG_kalman = zeros(grid_len_peaksAfterKalman_gt, 1);
    %sig_HR_kalman  = zeros(grid_len_peaksAfterKalman_gt, 1);

    idx_ecg = round(gt_peak_times .* fs_common)+1;
    idx_HR  = round(hr_peak_times .* fs_common)+1;
    %idx_HR_after_kalman = round(hr_peak_times_after_kalman .* fs_common)+1;

    sig_ECG(idx_ecg)=1;
    sig_HR(idx_HR)=1;

    sig_ECG_kalman(idx_ecg)=1;
    %sig_HR_kalman(idx_HR_after_kalman)=1;

    [res, lags_hr_gt] = xcorr(sig_HR, sig_ECG,2*fs_common);
    %[res_kalman, lags_hrkalman_gt] = xcorr(sig_HR_kalman, sig_ECG_kalman, 2*fs_common);
    
    [~,max_idx_normal]=max(res);
    %[~,max_idx_kalman]=max(res_kalman);


    lag_samples_normal = lags_hr_gt(max_idx_normal);
    %lag_samples_kalman = lags_hrkalman_gt(max_idx_kalman);

    delay_sec_normal = lag_samples_normal/fs_common;
    %delay_sec_kalman = lag_samples_kalman/fs_common;
    
    obj.mechanicalDelayMedian = delay_sec_normal;
   % obj.mechanicalDelayKalman = delay_sec_kalman;

    if (delay_sec_normal>0)
        sign=1;
    elseif (delay_sec_normal<0)
        sign=-1;
    end


end
function KalmanSmooth_BiDir(obj, Q_in, R_in, P_in)
    %KalmanSmooth_BiDir (Zero-Lag Smoother)
    %Runs the Kalman filter Forward, then flips the data and runs it Backward.
   % The average of both passes removes phase delay, massively improving correlation.

    %1. Parameter Handling
    if nargin < 2
        Q = obj.param_Q; R = obj.param_R; P = obj.param_P;
    else
        Q = Q_in; R = R_in; P = P_in;
    end

    %2. Safety Guards
    if isempty(obj.HrPeaks) || length(obj.HrPeaks) < 3
        obj.HrEstAfterKalman = [];
        return;
    end

    %3. Prepare Data (BPM Vector)
    raw_pk_times = unique(sort(obj.HrPeaks(:)), 'stable');
    ibi_meas = diff(raw_pk_times);
    hr_meas = 60 ./ ibi_meas; 

   % Filter out extreme outliers (40-200 BPM) before starting
    hr_meas(hr_meas < 40 | hr_meas > 200) = NaN;
    hr_meas = fillmissing(hr_meas, 'nearest'); % Fill gaps for smoother

   % --- PASS 1: FORWARD ---
    [x_fwd, ~] = obj.RunKalmanCore(hr_meas, Q, R, P);

    %--- PASS 2: BACKWARD ---
   % Flip the data to run backwards
    hr_meas_flipped = flipud(hr_meas);
    [x_bwd_flipped, ~] = obj.RunKalmanCore(hr_meas_flipped, Q, R, P);
    x_bwd = flipud(x_bwd_flipped); % Flip back to normal

    %--- COMBINE (Average) ---
    %This cancels the lag.
    hr_smooth = (x_fwd + x_bwd) / 2;

    obj.HrEstAfterKalman = hr_smooth;

   % --- Reconstruct Peaks (Optional Phase Locking) ---
    %Re-align peak times based on this smooth HR
    rec_peaks = zeros(length(hr_smooth)+1, 1);
    rec_peaks(1) = raw_pk_times(1);
    curr_time = raw_pk_times(1);
    for k = 1:length(hr_smooth)
        pred_ibi = 60 / hr_smooth(k);
        rec_peaks(k+1) = curr_time + pred_ibi;
        curr_time = rec_peaks(k+1);
    end

    %Simple alignment to nearest real peaks to prevent drift
   % (Only if real peaks are close)
    for k=1:length(rec_peaks)
        [min_val, min_idx] = min(abs(raw_pk_times - rec_peaks(k)));
        if min_val < 0.1 % If within 100ms
            rec_peaks(k) = raw_pk_times(min_idx);
        end
    end
    obj.HrPeaksAfterKalman = rec_peaks;
end
%Helper Function (Private or Public)
function [x_out, P_out] = RunKalmanCore(obj, z_in, Q, R, P)
%    The core mathematical loop, separated for reuse
    N = length(z_in);
    x_out = nan(N, 1);
    x = z_in(1); 
    if isnan(x), x = 70; end % Fallback start

    for k = 1:N
        z = z_in(k);

        %Predict
        x_pred = x;
        P_pred = P + Q;

       %Update (if valid)
        if isfinite(z)
            K = P_pred / (P_pred + R);
            x = x_pred + K * (z - x_pred);
            P = (1 - K) * P_pred;
        else
            x = x_pred;
            P = P_pred;
        end
        x_out(k) = x;
    end
    P_out = P;
end
%% PLOTS        
        % ---------------------------------------------------------
        % Plotting Functions
        % ---------------------------------------------------------
        
        %% Plot 1: Covariance Diagonal (Scatter)
       function h = plotCovDiag(obj,name,HrToCompare)
            h = []; % Default empty if no data
            % Check if estimates exist
            if isempty(obj.HrEst) || isempty(obj.HrGtEst)
                warning('HR Estimates are empty. Run FindRates() first.');
                return;
            end

            h = figure('Name', 'Covariance_Diag_Analysis', 'Color', 'w');
            
            % Sync lengths
            hold on;
           
            title(sprintf('Correlation Diag for GT vs %s - ID: %s, Scenario: %s',name, string(obj.ID), obj.sceneario));
            
            lengthmax= min(length(HrToCompare),length(obj.CorrGt));
            plot(obj.CorrGt(1:lengthmax),HrToCompare(1:lengthmax),"DisplayName",' Median Ground truth to RADAR heart rate'...
                ,'Marker','*','LineStyle','none')
            xlabel('Ground Truth HR (Post Median Filter)');
            ylabel(['RADAR Estimated HR of ', name]);
            
            plot((30:1:150),(30:1:150));
            rho = corr(obj.CorrGt(:), HrToCompare(:), 'Rows', 'complete');
            legend(sprintf('Measurements (r = %.2f)', rho), 'Location', 'best');
        end       
       
       %% Plot 2: HR Peaks and Signals
       function h = plotHrpeaks(obj)
            h = figure('Name', 'Radar_ECG_Peaks', 'Color', 'w');
            
            % Subplot 1: Radar Signal
            ax(1) = subplot(2,1,1);
            hold on;
            title(sprintf('Radar HR Filters & Peak Finder - ID: %s, Scenario: %s', string(obj.ID), obj.sceneario));
            
            plot(obj.vTimeNew, obj.HrSignal*1e4, 'b', 'DisplayName', 'Filtered Radar');
            
            % Interpolate to find amplitude at peak times
            if ~isempty(obj.HrPeaks)
                peakAmps = interp1(obj.vTimeNew, obj.HrSignal, obj.HrPeaks); 
                peakAmpsAfterKalman = interp1(obj.vTimeNew, obj.HrSignal, obj.HrPeaksAfterKalman); 
                plot(obj.HrPeaks, peakAmps*1e4, 'r*', 'MarkerSize', 8, 'DisplayName', 'Radar Peaks'); 
                plot(obj.HrPeaksAfterKalman, peakAmpsAfterKalman*1e4, 'go', 'MarkerSize', 8, 'DisplayName', 'Radar Peaks after kalman filter'); 
            end
            
            ylabel('Amp (scaled)');
            xlabel('Time (s)'); 
            legend('show', 'Location', 'best'); grid on; hold off;

            % Subplot 2: ECG Reference
            ax(2) = subplot(2,1,2);
            hold on;
            title(sprintf('ECG Reference Signal - ID: %s, Scenario: %s', string(obj.ID), obj.sceneario));
            
            plot(obj.vTimeOriginal, obj.ecg_gt, 'k', 'DisplayName', 'ECG Signal');
            
            if ~isempty(obj.ecgPeaks)
                peakAmpsEcg = interp1(obj.vTimeOriginal, obj.ecg_gt, obj.ecgPeaks);
                plot(obj.ecgPeaks, peakAmpsEcg, 'r*', 'MarkerSize', 8, 'DisplayName', 'ECG Peaks');
            end
            
            xlabel('Time (s)'); ylabel('Amp');
            legend('show', 'Location', 'best'); grid on; hold off;
            
            linkaxes(ax, 'x');
        end

        %% Plot 3: Bland-Altman
       function h = plotBA(obj, name,HrToCompare)
            h = []; % Default empty if no data
            % Check if estimates exist
            if isempty(obj.HrEst) || isempty(obj.HrGtEst)
                warning('HR Estimates are empty. Run FindRates() first.');
                return;
            end

            h = figure('Name', 'Bland_Altman_Analysis', 'Color', 'w');
            if ~isempty(obj.CorrGt) && ~isempty(HrToCompare)
                if exist('BlandAltman', 'file')
                   maxlen=min(length(HrToCompare),length(obj.CorrGt));
                   BlandAltman(HrToCompare(1:maxlen),obj.CorrGt(1:maxlen), 2, 0);
                end
                title(sprintf('Bland-Altman of %s - ID: %s, Scenario: %s',name,string(obj.ID), obj.sceneario));
            end
        end
       
       %% Plot 4: Respiration Rates (Trends)
       function h = plotRrRates(obj)
            h = figure('Name', 'Respiration_Rates', 'Color', 'w');
            hold on;
            title(sprintf('Respiration Rate Comparison - ID: %s, Scenario: %s', string(obj.ID), obj.sceneario));
           
            if ~isempty(obj.RrEst)
                plot(obj.RrEst, 'b.-', 'LineWidth', 1.5, 'DisplayName', 'Radar Respiration');
            end

            if ~isempty(obj.RrGtEst)
                 plot(obj.RrGtEst, 'r.-', 'LineWidth', 1.5, 'DisplayName', 'TFM Respiration');           
            end
            ylabel('Breaths Per Minute'); 
            xlabel('Window Index / Time');
            legend('show', 'Location', 'best'); 
            grid on; hold off;
        end

        %% Plot 5: Respiration Signals (Time Domain)
       function h = plotRrSignals(obj)
            h = figure('Name', 'Respiration_Signals_Comparison', 'Color', 'w');
            
            ax = [];
            
            % Subplot 1: Radar Respiration Signal
            ax(1) = subplot(2,1,1);
            hold on;
            title(sprintf('Radar Respiration Signal - ID: %s, Scenario: %s', string(obj.ID), obj.sceneario));
            
            plot(obj.vTimeNew, obj.RrSignal, 'b-', 'DisplayName', 'Respiration Signal');
            if ~isempty(obj.RrPeaks)
                peakAmps = interp1(obj.vTimeNew, obj.RrSignal, obj.RrPeaks);
                plot(obj.RrPeaks, peakAmps, 'k*', 'MarkerSize', 8, 'DisplayName', 'Radar Peaks');
            end
            ylabel('Rel. Distance(mm)');
            xlabel('Time (s)'); % ADDED
            legend('show', 'Location', 'best'); grid on; hold off;

            % Subplot 2: GT Respiration Signal
            ax(2) = subplot(2,1,2);
            hold on;
            title(sprintf('GT Respiration Signal (TFM) - ID: %s, Scenario: %s', string(obj.ID), obj.sceneario));

            % Create time vector for GT if not stored
            time_respiration = 1/100:1/100:length(obj.resp_gt)/100;
            
            plot(time_respiration, obj.resp_gt, 'r-', 'DisplayName', 'GT Signal');
            if ~isempty(obj.Rrpeaks_gt)
                peakAmps = interp1(time_respiration, obj.resp_gt, obj.Rrpeaks_gt);
                plot(obj.Rrpeaks_gt, peakAmps , 'k*', 'MarkerSize', 8, 'DisplayName', 'GT Peaks');
            end
            
            ylabel('Rel. Distance(mm)');
            xlabel('Time (s)');
            legend('show', 'Location', 'best'); grid on; hold off;          
            
            linkaxes(ax, 'x');
        end    
       
       %% Plot 6: Dashboard (Summary)
       function h = plotDashBoard(obj,name,HrToCompare,peaksToCompare,timeToCompare) 
            h = figure('Name', 'Summary_Dashboard', 'Units', 'normalized', 'Position', [0.1 0.1 0.6 0.8], 'Color', 'w');
            
            ax_link = [];

            % 1. Radar Heart signal WITH PEAKS
            ax_link(1) = subplot(4,1,1);
            hold on;
            plot(obj.vTimeNew, obj.HrSignal, 'b', 'DisplayName', 'Radar Signal');
            if ~isempty(obj.HrPeaks)
                peakAmps = interp1(obj.vTimeNew, obj.HrSignal, obj.HrPeaks); 
                plot(obj.HrPeaks, peakAmps, 'r*', 'MarkerSize', 8, 'DisplayName', 'Radar Peaks');
                peakAmpsToCompare = interp1(obj.vTimeNew, obj.HrSignal, peaksToCompare);
                 plot(peaksToCompare, peakAmpsToCompare, 'go', 'MarkerSize', 8, 'DisplayName', 'Radar Peaks');
            end
            title(sprintf('Radar Signal - ID: %s, Scenario: %s', string(obj.ID), obj.sceneario));
            ylabel('Amp'); 
            xlabel('Time (s)'); % ADDED
            legend('show', 'Location', 'best'); grid on; axis tight; hold off;

            % 2. ECG reference WITH PEAKS
            ax_link(2) = subplot(4,1,2);
            hold on;
            plot(obj.vTimeOriginal, obj.ecg_gt, 'k', 'DisplayName', 'ECG Signal');
            if ~isempty(obj.ecgPeaks)
                peakAmpsEcg = interp1(obj.vTimeOriginal, obj.ecg_gt, obj.ecgPeaks);
                plot(obj.ecgPeaks, peakAmpsEcg, 'r*', 'MarkerSize', 8, 'DisplayName', 'QRS Peaks');
            end
            title(sprintf('ECG Reference - ID: %s, Scenario: %s', string(obj.ID), obj.sceneario));
            ylabel('Amp'); 
            xlabel('Time (s)'); % ADDED
            legend('show', 'Location', 'best'); grid on; axis tight; hold off;

            % 3. Heart Rate Comparison (BPM vs Time)
            ax_link(3) = subplot(4,1,3);
            hold on;
            if ~isempty(obj.HrGtEst)
                time_gt_bpm = obj.ecgPeaks(2:end);
                plot(time_gt_bpm, obj.HrGtEstAfterMedian, 'r.-', 'LineWidth', 1.5, 'DisplayName', 'ECG GT (after Median)');
            end
            
            plot(timeToCompare,HrToCompare, 'b.--', 'LineWidth', 1.2, 'DisplayName', name);
            title(sprintf('Heart Rate (BPM) - ID: %s, Scenario: %s', string(obj.ID), obj.sceneario));
            xlabel('Time (s)'); ylabel('BPM'); legend('show', 'Location', 'best'); grid on; axis tight; hold off;

            % 4. Bland-Altman (NOT LINKED)
            subplot(4,1,4);
            if ~isempty(obj.CorrGt) && ~isempty(HrToCompare)
                if exist('BlandAltman', 'file')
                   maxlen=min(length(HrToCompare),length(obj.CorrGt));
                   BlandAltman(HrToCompare(1:maxlen),obj.CorrGt(1:maxlen), 2, 0);
                end
                title(sprintf('Bland-Altman - ID: %s, Scenario: %s', string(obj.ID), obj.sceneario));
            end

            linkaxes(ax_link, 'x');            
        end
       
       %% Plot 7: Examples (Comparison)
       function h = plot_examples(obj)
         h = gobjects(0); % Initialize graphics array

         % --- FIGURE 1: Peaks Comparison on HrSignal & ECG ---
         h(end+1) = figure('Name', 'Peaks_Comparison', 'Color', 'w');
         
         ax_link = [];

         % Subplot 1
         ax_link(1) = subplot(2,1,1);
         hold on;
         title(sprintf('Radar Signal: Initial vs Final Peaks - ID: %s, Scenario: %s', string(obj.ID), obj.sceneario));
         plot(obj.vTimeNew, obj.HrSignal, 'b', 'DisplayName', 'Radar Signal');
         
         if ~isempty(obj.HrPeaks)
            amps = interp1(obj.vTimeNew, obj.HrSignal, obj.HrPeaks);
            plot(obj.HrPeaks, amps, 'g*', 'MarkerSize', 6, 'DisplayName', 'Initial Peaks');
         end
         if ~isempty(obj.HrPeaksAfterKalman)
            amps = interp1(obj.vTimeNew, obj.HrSignal, obj.HrPeaksAfterKalman);
            plot(obj.HrPeaksAfterKalman, amps, 'ro', 'MarkerSize', 8, 'DisplayName', 'After kalman Peaks');
         end
         ylabel('Amplitude'); % ADDED
         xlabel('Time (s)');  % ADDED
         legend('show', 'Location', 'best'); grid on; hold off;

         % Subplot 2
         ax_link(2) = subplot(2,1,2);
         hold on;
         title('ECG Signal: GT Peaks vs Correlated Radar Peaks');
         plot(obj.vTimeOriginal, obj.ecg_gt, 'k', 'DisplayName', 'ECG Signal');
         
         if ~isempty(obj.ecgPeaks)
             amps = interp1(obj.vTimeOriginal, obj.ecg_gt, obj.ecgPeaks);
             plot(obj.ecgPeaks, amps, 'bo', 'MarkerSize', 6, 'DisplayName', 'ECG Peaks (GT)');
         end       
         ylabel('Amplitude'); % ADDED
         xlabel('Time (s)');  % ADDED
         legend('show', 'Location', 'best'); grid on; hold off;
         linkaxes(ax_link, 'x');

         % --- FIGURE 2: Heart Rate Comparisons ---
         h(end+1) = figure('Name', 'Heart_Rate_Comparison', 'Color', 'w');
         
         ax_hr = [];

         % Subplot 1
         ax_hr(1) = subplot(3,1,1); hold on; grid on;
         title(sprintf('HR Comparison: Radar Raw vs GT (Median Fit) - ID: %s, Scenario: %s', string(obj.ID), obj.sceneario));
         plot(obj.HrEst, 'b', 'DisplayName', 'Radar Est - no processing');
         plot(obj.CorrGt, 'r--','DisplayName', 'GT(Median&fit)');
         ylabel('Heart Rate (BPM)'); % ADDED
         xlabel('Beat Index (Approx)'); % ADDED
         legend('show', 'Location', 'best');

         % Subplot 2
         ax_hr(2) = subplot(3,1,2); hold on; grid on;
         title('HR after Kalman (time fitted) vs GT (time fitted)');
          if ~isempty(obj.CorrKalmanHr_on_gt_time)
             plot(obj.CorrKalmanHr_on_gt_time, 'b','DisplayName', 'Radar Est after Kalman(time fitted)');
          end
         plot(obj.CorrGt, 'r--',  'DisplayName', 'GT(Median&fit)');
         ylabel('Heart Rate (BPM)'); % ADDED
         xlabel('Time (s) [GT Grid]'); % ADDED
         legend('show', 'Location', 'best');

         %Subplot 3
         ax_hr(3) = subplot(3,1,3); hold on; grid on;
         title('HR after Median vs GT (time fitted)');
         if ~isempty(obj.CorrMedianHr_on_gt_time)
             plot(obj.CorrMedianHr_on_gt_time, 'b','DisplayName', 'Radar Est After Median');
         end
         plot(obj.CorrGt, 'r--','DisplayName', 'GT');
         ylabel('Heart Rate (BPM)'); % ADDED
         xlabel('Time (s) [GT Grid]'); % ADDED
         legend('show', 'Location', 'best');

         linkaxes(ax_hr, 'xy');
     end

       %% Plot 8: Error Analysis
       function h = plotErrors(obj)
            h = figure('Name', 'HR_Error_Analysis', 'Color', 'w');
            
            % --- Subplot 1: Raw Error vs HR ---
            subplot(3,1,1);
            hold on; grid on;
            if ~isempty(obj.rawErr2Hr)
                plot(obj.rawErr2Hr(:,1), obj.rawErr2Hr(:,2), 'ko', 'MarkerFaceColor', 'b', 'MarkerSize', 4);
            end
            yline(0, 'r--', 'LineWidth', 1.5);
            % Updated Title with ID and Scenario
            title(sprintf('Raw Error (Bias) vs ECG HR - ID: %s, Scenario: %s', string(obj.ID), obj.sceneario));
            ylabel('Error');
            xlabel('ECG Heart Rate (BPM)');
            legend('Raw Error', 'Zero Line', 'Location', 'best');
            
            % --- Subplot 2: Absolute Error vs HR (The "MAE" request) ---
            subplot(3,1,2);
            hold on; grid on;
            if ~isempty(obj.mae2Hr)
                plot(obj.mae2Hr(:,1), obj.mae2Hr(:,2), 'ko', 'MarkerFaceColor', 'g', 'MarkerSize', 4);
            end
            % Updated Title with ID and Scenario
            title(sprintf('Absolute Error (Magnitude) vs ECG HR - ID: %s, Scenario: %s', string(obj.ID), obj.sceneario));
            ylabel('|Error| ');
            xlabel('ECG Heart Rate (BPM)');
            yline(mean(obj.mae2Hr(:,2)), 'm--', 'DisplayName', 'Mean (MAE)');
            legend('Abs Error', 'MAE Level', 'Location', 'best');

            % --- Subplot 3: Squared Error vs HR (MSE) ---
            subplot(3,1,3);
            hold on; grid on;
            if ~isempty(obj.mse2Hr)
                plot(obj.mse2Hr(:,1), obj.mse2Hr(:,2), 'ko', 'MarkerFaceColor', 'r', 'MarkerSize', 4);
            end
            % Updated Title with ID and Scenario
            title(sprintf('Squared Error (Outliers) vs ECG HR - ID: %s, Scenario: %s', string(obj.ID), obj.sceneario));
            xlabel('ECG Heart Rate (BPM)');
            ylabel('Error^2');
            
            % Print Stats
            fprintf('------------------------------------------------\n');
            fprintf('Error Analysis for ID: %s, Scenario: %s\n', obj.ID, obj.sceneario);
            fprintf('Mean Absolute Error (MAE) -> RAW: %.2f | Fitted: %.2f\n', obj.maeRaw, obj.maeFitted);
            fprintf('Mean Squared Error (MSE)  -> RAW: %.2f | Fitted: %.2f\n', obj.mseRaw, obj.mseFitted);
            fprintf('------------------------------------------------\n');
       end
  %% NEW FUNCTION
  function PlotHrCovAndBA(obj)
% PlotHrCovAndBA
% Creates two figures:
%   (1) 2x2 subplots: "CovDiag" style plots for:
%       - HrEst vs HrGtEst
%       - HrEstAfterMedian vs HrGtEst
%       - HrEstAfterKalman vs HrGtEst
%   (2) 2x2 subplots: Bland-Altman plots for the same three comparisons
%
% Notes:
%   - GT is obj.HrGtEst (as requested)
%   - Vectors are cut to the same length and NaNs are removed consistently.

    gt = obj.HrGtEst;

    sigs = {obj.HrEst, obj.HrEstAfterMedian, obj.HrEstAfterKalman};
    names = {'Raw HR vs GT', 'Median HR vs GT', 'Kalman HR vs GT'};

    % ===================== FIGURE 1: CovDiag-style =====================
    figure('Name', sprintf('CovDiag HR vs GT | ID %s | %s', string(obj.ID), string(obj.sceneario)));

    for k = 1:3
        subplot(2,2,k);

        est = sigs{k};
        a = est(:); b = gt(:);

        n = min(numel(a), numel(b));
        a = a(1:n); b = b(1:n);

        v = isfinite(a) & isfinite(b) & ~isnan(a) & ~isnan(b);
        a = a(v); b = b(v);

        if numel(a) < 8
            axis off;
            title([names{k} ' (insufficient data)']);
            continue;
        end

        % Scatter
        plot(b, a, '.', 'MarkerSize', 8); hold on; grid on;
        xlabel('GT HR'); ylabel('Est HR');
        title(names{k});

        % y=x reference line
        mn = min([a; b]);
        mx = max([a; b]);
        plot([mn mx], [mn mx], 'k--', 'LineWidth', 1);

        % --- Covariance ellipse (1-sigma) in (GT,Est) space ---
        % Data matrix: columns [GT, Est]
        X = [b, a];
        mu = mean(X, 1);
        C  = cov(X, 1); % population cov

        % Eigen-decomp
        [V,D] = eig(C);
        d = diag(D);
        [d,idx] = sort(d, 'descend');
        V = V(:,idx);

        % 1-sigma ellipse
        t = linspace(0, 2*pi, 200);
        circ = [cos(t); sin(t)];
        A = V * diag(sqrt(max(d,0)));   % sqrt eigenvalues
        ell = (A * circ).';
        ell(:,1) = ell(:,1) + mu(1);
        ell(:,2) = ell(:,2) + mu(2);

        plot(ell(:,1), ell(:,2), 'LineWidth', 1.5);

        % --- Basic stats annotation ---
        r = corr(a, b);
        diffv = a - b;
        bias = mean(diffv);
        rmse = sqrt(mean(diffv.^2));

        txt = sprintf('N=%d  r=%.3f\nbias=%.2f  rmse=%.2f', numel(a), r, bias, rmse);
        xlim([mn mx]); ylim([mn mx]);
        text(mn + 0.02*(mx-mn), mx - 0.10*(mx-mn), txt, 'FontSize', 9, 'BackgroundColor', 'w');
        hold off;
    end

    subplot(2,2,4);
    axis off;
    text(0,0.85, sprintf('ID: %s', string(obj.ID)), 'FontWeight','bold');
    text(0,0.65, sprintf('Scenario: %s', string(obj.sceneario)));
    text(0,0.45, 'CovDiag-style: scatter + y=x + 1σ covariance ellipse');
    text(0,0.25, 'All compared to GT = HrGtEst');
    text(0,0.05, 'Vectors trimmed to same length; NaNs removed');

    % ===================== FIGURE 2: Bland-Altman =====================
    figure('Name', sprintf('Bland-Altman HR vs GT | ID %s | %s', string(obj.ID), string(obj.sceneario)));

    for k = 1:3
        subplot(2,2,k);

        est = sigs{k};
        a = est(:); b = gt(:);

        n = min(numel(a), numel(b));
        a = a(1:n); b = b(1:n);

        v = isfinite(a) & isfinite(b) & ~isnan(a) & ~isnan(b);
        a = a(v); b = b(v);

        if numel(a) < 8
            axis off;
            title([names{k} ' (insufficient data)']);
            continue;
        end

        m  = 0.5*(a + b);      % mean
        d  = (a - b);          % difference (Est - GT)
        md = mean(d);
        sd = std(d, 0);

        loa1 = md - 1.96*sd;
        loa2 = md + 1.96*sd;

        plot(m, d, '.', 'MarkerSize', 8); hold on; grid on;
        yline(md,  'k-',  'LineWidth', 1.5);
        yline(loa1,'k--', 'LineWidth', 1.0);
        yline(loa2,'k--', 'LineWidth', 1.0);

        xlabel('Mean (Est, GT)'); ylabel('Est - GT');
        title(names{k});

        txt = sprintf('N=%d\nmean=%.2f\nLoA=[%.2f, %.2f]', numel(a), md, loa1, loa2);
        xm = min(m); xM = max(m);
        ym = min(d); yM = max(d);
        text(xm + 0.02*(xM-xm), yM - 0.10*(yM-ym), txt, 'FontSize', 9, 'BackgroundColor', 'w');

        hold off;
    end

    subplot(2,2,4);
    axis off;
    text(0,0.85, sprintf('ID: %s', string(obj.ID)), 'FontWeight','bold');
    text(0,0.65, sprintf('Scenario: %s', string(obj.sceneario)));
    text(0,0.45, 'Bland-Altman: (Est-GT) vs mean, with mean & ±1.96σ');
    text(0,0.25, 'All compared to GT = HrGtEst');
    text(0,0.05, 'Vectors trimmed to same length; NaNs removed');

end




        %% Save Figures
       function [] = saveFigures(obj, figHandles, saveDir)
             % Set default directory if not provided
            if nargin < 3
                saveDir = 'SavedAnalysisFigures';
            end

            % Filter out invalid handles
            if isempty(figHandles)
                fprintf('No figure handles provided for ID %s.\n', string(obj.ID));
                return;
            end

            % Create directory if it doesn't exist
            if ~exist(saveDir, 'dir')
                mkdir(saveDir);
                fprintf('Created save directory: %s\n', saveDir);
            end

            fprintf('Saving %d figures to: %s\n', length(figHandles), saveDir);

            for i = 1:length(figHandles)
                hFig = figHandles(i);
                
                % Use isgraphics to ensure it's a valid figure
                if ~isgraphics(hFig)
                    continue; 
                end

                % Get Name
                figName = get(hFig, 'Name');
                if isempty(figName)
                    figName = sprintf('Figure%d', get(hFig, 'Number'));
                end

                % Clean filename
                safeFigName = strrep(figName, ' ', '_');
                safeFigName = regexprep(safeFigName, '[^a-zA-Z0-9_]', '');
                
                % Construct full filename: ID_Scenario_FigName
                baseFileName = sprintf('%s_%s_%s', string(obj.ID), obj.sceneario, safeFigName);
                
                %pngFileName = fullfile(saveDir, [baseFileName, '.png']);
                figFileName = fullfile(saveDir, [baseFileName, '.fig']);

                try
                    %saveas(hFig, pngFileName, 'png'); 
                    saveas(hFig, figFileName, 'fig');
                    fprintf('  -> Saved: %s\n', baseFileName);
                catch ME
                    fprintf(2, 'Warning: Could not save %s. Error: %s\n', baseFileName, ME.message);
                end
             end
        end
         %% Plot All with Selection Flags
       function [] = PlotAll(obj, bsave, saveDir,name,HrToCompare,peaksToCompare,timeToCompare, options) 
            arguments
                obj
                bsave (1,1) logical = false
                saveDir (1,1) string = 'SavedAnalysisFigures'
                name {mustBeText} = 'INPUT NAME TO PlotALL' %NEW, REORDER CALLS TO PLOTALL
                HrToCompare = obj.HrEst ;
                peaksToCompare = obj.HrPeaks ;
                timeToCompare = obj.vTimeNew ;
                options.plot_HrPeaks (1,1) logical = true
                options.plot_RrRates (1,1) logical = true
                options.plot_RrSignals (1,1) logical = true
                options.plot_BA (1,1) logical = true
                options.plot_DashBoard (1,1) logical = true
                options.plot_Errors (1,1) logical = true
                options.plot_CorrAxis (1,1) logical = true
                options.plot_examples (1,1) logical = true 
            end
            figHandles = gobjects(0); % Initialize empty graphics array

            % 1. HR Peaks
            if options.plot_HrPeaks
                h1 = obj.plotHrpeaks();
                if isgraphics(h1), figHandles(end+1) = h1; end
            end
            
            % 2. Respiration Rates
            if options.plot_RrRates
                h2 = obj.plotRrRates();
                if isgraphics(h2), figHandles(end+1) = h2; end
            end

            % 3. Respiration Signals
            if options.plot_RrSignals
                h3 = obj.plotRrSignals();
                if isgraphics(h3), figHandles(end+1) = h3; end
            end

            % 4. Bland Altman (only if data exists)
            if options.plot_BA
                h4 = obj.plotBA(name,HrToCompare);
                if isgraphics(h4), figHandles(end+1) = h4; end
            end
             if options.plot_CorrAxis %TODO is this right?
                h41 = obj.plotCovDiag(name,HrToCompare);
                if isgraphics(h41), figHandles(end+1) = h41; end
            end

            % 5. Dashboard
            if options.plot_DashBoard
                h5 = obj.plotDashBoard(name,HrToCompare,peaksToCompare,timeToCompare);
                if isgraphics(h5), figHandles(end+1) = h5; end
            end

            % 6. Error Analysis
            if options.plot_Errors
                h6 = obj.plotErrors();
                if isgraphics(h6), figHandles(end+1) = h6; end
            end

            % 7. Examples (New Figures)
            if options.plot_examples
                h7 = obj.plot_examples();
                if ~isempty(h7)
                     for k=1:length(h7)
                         if isgraphics(h7(k))
                            figHandles(end+1) = h7(k); 
                         end
                     end
                end
            end

            % Save if requested
            if bsave
                obj.saveFigures(figHandles, saveDir);
                close all; 
            end
        end
    end
    methods (Access = private)
        function y_interp = replaceOutliersByInterp(obj, y, outlierIdxOrMask) %#ok<INUSL>
            % Replace outliers by interpolation over index/time.
            % y: vector
            % outlierIdxOrMask: logical mask same length as y OR index vector
    
            y = y(:);
            N = numel(y);
            t = (1:N).';
    
            % Build logical mask
            if islogical(outlierIdxOrMask)
                bad = outlierIdxOrMask(:);
                if numel(bad) ~= N
                    error('replaceOutliersByInterp:MaskSize', ...
                          'Logical mask must be same length as y.');
                end
            else
                bad = false(N,1);
                idx = outlierIdxOrMask(:);
                idx = idx(idx>=1 & idx<=N & isfinite(idx));
                bad(idx) = true;
            end
    
            good = ~bad;
    
            % Need at least 2 good points to interpolate
            if nnz(good) < 2
                y_interp = y;  % nothing we can do robustly
                return;
            end
    
            y_interp = y;
    
            % Interpolate only the bad samples
            y_interp(bad) = interp1(t(good), y(good), t(bad), 'pchip', 'extrap');
        end
    end

end
