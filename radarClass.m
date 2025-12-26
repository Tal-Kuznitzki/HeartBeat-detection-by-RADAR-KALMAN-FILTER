classdef radarClass < handle
    properties
        %instance details
        ID
        sceneario
        %basic signals
        radar_i
        radar_q
        radar_dist
        ecg_gt
        resp_gt
        
        %proccessed signals
        radar_decimated
        HrSignal
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
        correlation_misses
        correlation_excess


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
                scenario (1,1) string {mustBeMember(scenario, ["Resting","Valsalva","Apnea","Tilt-down","Tilt-up"])} 
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
            sprintf('Set new radar fs to %d',fs)
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
            obj.RrPeaks = obj.RrPeaks / obj.fs_new;
            obj.ecgPeaks = obj.ecgPeaks / obj.fs_ecg;
            obj.Rrpeaks_gt = obj.Rrpeaks_gt /(100);
        end

        % finding the rates
        function FindRates(obj)
            obj.HrEst = 60 ./  diff(obj.HrPeaks);
            obj.HrGtEst = 60 ./  diff(obj.ecgPeaks); 
            obj.RrEst = 60 ./  diff(obj.RrPeaks);
            obj.RrGtEst = 60 ./ diff(obj.Rrpeaks_gt);            
        end
        function MedianHr(obj,size)
            if(nargin<2)
                size=10;
            end
          gtHr=obj.HrGtEst(:);
          calculatedHr=obj.HrEst(:);
          minlen=min(length(gtHr),length(calculatedHr));
          gtHr=obj.HrGtEst(1:minlen);
          calculatedHr=obj.HrEst(1:minlen);
          gt_after_median_filt=medfilt1(gtHr,size);
          calculatedHr_after_median_filt=medfilt1(calculatedHr,size);


          %FOR RADAR SIGNAL
         indx = find(calculatedHr > 1.2 * calculatedHr_after_median_filt | calculatedHr < 0.8 * calculatedHr_after_median_filt); 
            
          % hr_replaced_with_median_in_outliers = calculatedHr;
          % hr_replaced_with_median_in_outliers(indx) =  calculatedHr_after_median_filt(indx);
          hr_replaced_with_median_in_outliers = obj.replaceOutliersByInterp(calculatedHr, indx);

          indx = find(gtHr > 1.2 * gt_after_median_filt | gtHr < 0.8 * gt_after_median_filt); 
            
          % gt_replaced_with_median_in_outliers = gtHr;
          % gt_replaced_with_median_in_outliers(indx) =  gt_after_median_filt(indx);
          gt_replaced_with_median_in_outliers = obj.replaceOutliersByInterp(gtHr, indx);

          
        hr_replaced_with_median_in_outliers=hr_replaced_with_median_in_outliers(:);
        gt_replaced_with_median_in_outliers=gt_replaced_with_median_in_outliers(:);
       % cor_after_median_Filter_on_both_signals = corr(hr_replaced_with_median_in_outliers,gt_replaced_with_median_in_outliers)
        obj.HrEstAfterMedian = hr_replaced_with_median_in_outliers;
        obj.HrGtEstAfterMedian = gt_replaced_with_median_in_outliers;

        % correlation_value = R(1,2);
        % fprintf('The correlation coefficient is: %.4f\n', correlation_value);
        % 
        % 
        % 
        % 
        % 
        % 
        %   h = figure('Name', 'diff', 'Color', 'w');
        %   legend('show', 'Location', 'best');
        %   subplot(3,1,1);
        %    plot(gtHr,'Color', 'b', 'DisplayName', 'GT');
        %    hold on;
        %    plot(gt_after_median_filt,'Color', 'r', 'DisplayName', 'gt filtered');
        %    plot(gt_replaced_with_median_in_outliers,'Color', 'g', 'DisplayName', 'after median fix');     
        %    legend('show', 'Location', 'best');
        %    subplot(3,1,2);
        %    plot(calculatedHr,'Color', 'b', 'DisplayName', 'Radar Signal');
        %    hold on;
        %    plot(calculatedHr_after_median_filt,'Color', 'r', 'DisplayName', 'radar filtered');
        %    plot(hr_replaced_with_median_in_outliers,'Color', 'g', 'DisplayName', 'after median fix');
        %    legend('show', 'Location', 'best');
        %    subplot(3,1,3);
        %    plot(gt_replaced_with_median_in_outliers,'Color', 'r', 'DisplayName', 'gt post fix');
        %    hold on;
        %    plot(hr_replaced_with_median_in_outliers,'Color', 'g', 'DisplayName', 'radar post fix');
        %    legend('show', 'Location', 'best');
        % 
        % 
        %   linkaxes
        end
        function [] = CalcError(obj)
            min_len = min(length(obj.HrGtEst), length(obj.HrEstAfterKalman));
            v1=obj.HrEstAfterKalman(1:min_len);
            v2=obj.HrGtEst(1:min_len);
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

%% KALMAN 

%OLD KALMAN IMPL

% function KalmanFilterBeats(obj, b_both)
% % KalmanFilterBeats
% % Applies a Kalman filter on beat-to-beat HR derived from obj.HrPeaks and:
% %   1) saves corrected peak times in obj.HrPeaksAfterKalman
% %   2) saves HR estimate (bpm) in obj.HrEstAfterKalman (beat-rate sequence)
% % If b_both==1: also applies a KF directly to obj.HrEst (if exists) and
% % uses that filtered HR instead (or in addition) for HrEstAfterKalman.
% %
% % Notes:
% % - This does NOT filter the waveform; it filters the HR state.
% % - Peak reconstruction uses: peaks_rec = t0 + [0; cumsum(60./HRhat)]
% %   where t0 is the first peak time.
% 
% if nargin < 2 || isempty(b_both)
%     b_both = 0;
% end
% 
% % ---- guards ----
% if isempty(obj.HrPeaks)
%     obj.HrPeaksAfterKalman = [];
%     obj.HrEstAfterKalman   = [];
%     return;
% end
% 
% fs = obj.fs_new;
% 
% % Ensure column, sorted, unique
% pk = unique(obj.HrPeaks(:), 'stable'); %stable keeps the order 
% pk = sort(pk);
% 
% if numel(pk) < 2
%     obj.HrPeaksAfterKalman = pk;
%     obj.HrEstAfterKalman   = [];
%     return;
% end
% 
% t = pk ;                 % seconds
% ibi = diff(t);               % seconds
% hr_meas = 60 ./ ibi;         % bpm (length N-1)
% 
% % ---- basic sanity (optional but recommended) ----
% HRmin = 40; HRmax = 200;
% good = isfinite(hr_meas) & hr_meas >= HRmin & hr_meas <= HRmax;
% % keep alignment: we'll KF-update only where good, but still predict through all beats
% 
% % ---- Kalman init (state: [HR; HRdot]) ----
% % Initial HR from median of first few valid beats
% hr0 = hr_meas(good);
% if isempty(hr0)
%     HR0 = 75;
%     disp("ERROR");
% else
%     HR0 = median(hr0(1:min(5,numel(hr0))));
%     HR0 = median(hr0);
% end
% HR0 = min(max(HR0, HRmin), HRmax);
% HR0_d=0;
% x = [HR0; HR0_d];                % [HR; HRdot], bpm and bpm/s
% P = diag([25, 9]);
% 
% % Tunables (these are decent starting points)
% Q = diag([0.8^2, 0.6^2]);     % process noise
% R0 = 4.0^2;                  % measurement noise (bpm^2)
% gateSig = 3.0;               % innovation gate in sigma
% 
% H = [1 0];
% 
% hr_hat = nan(size(hr_meas));
% acc    = false(size(hr_meas));
% 
% for k = 1:numel(hr_meas)
%     dt = ibi(k);
%     if ~isfinite(dt) || dt <= 0
%         dt = 0.8; % fallback
%     end
% 
%     % Predict
%     F = [1 dt; 0 1];
%     x = F*x;
%     P = F*P*F.' + Q;
% 
%     if good(k)
%         z = hr_meas(k);
% 
%         % Innovation + gate
%         y = z - H*x;
%         S = H*P*H.' + R0;
% 
%         if (y*y) <= (gateSig^2)*S
%             K = (P*H.') / S;
%             x = x + K*y;
%             P = (eye(2) - K*H) * P;
%             acc(k) = true;
%         end
%     end
% 
%     % Clamp HR estimate
%     x(1) = min(max(x(1), HRmin), HRmax);
%     hr_hat(k) = x(1);
% end
% 
% % ---- Optionally KF-filter an existing HR time-series (b_both == 1) ----
% % If obj.HrEst exists and matches length, filter it too and prefer it.
% if b_both == 1 && isprop(obj,'HrEst') && ~isempty(obj.HrEst)
%     hr_in = obj.HrEst(:);
% 
%     % If hr_in length matches hr_hat length, we can filter it beat-wise.
%     % Otherwise we just ignore it here (you can adapt later by resampling).
%     if numel(hr_in) == numel(hr_hat)
%         % Simple 1D KF on HR only (random walk)
%         x1 = hr_in(1);
%         P1 = 25;
%         Q1 = 0.8^2;
%         R1 = 4.0^2;
% 
%         hr1 = nan(size(hr_in));
%         for k = 1:numel(hr_in)
%             % predict
%             x1 = x1;
%             P1 = P1 + Q1;
% 
%             z = hr_in(k);
%             if isfinite(z) && z >= HRmin && z <= HRmax
%                 % update
%                 K1 = P1 / (P1 + R1);
%                 x1 = x1 + K1*(z - x1);
%                 P1 = (1 - K1)*P1;
%             end
%             x1 = min(max(x1, HRmin), HRmax);
%             hr1(k) = x1;
%         end
% 
%         % Prefer filtered HR estimate (replace hr_hat)
%         hr_hat = hr1;
%     end
% end
% 
% % ---- Save HR estimate after Kalman ----
% % hr_hat is beat-to-beat (between peaks). Save as-is:
% obj.HrEstAfterKalman = hr_hat;
% 
% % ---- Reconstruct peaks from filtered HR (requires first peak time) ----
% t0 = t(1);                          % seconds, timestamp of first beat
% ibi_hat = 60 ./ hr_hat;             % seconds per beat
% 
% % Handle any NaNs (no update segments): fill with nearest valid
% if any(~isfinite(ibi_hat))
%     ibi_hat = fillmissing(ibi_hat, 'nearest');
% end
% 
% t_rec = t0 + [0; cumsum(ibi_hat)];
% pk_rec = round(t_rec * fs);
% 
% % Enforce monotonic increasing and minimum distance of 1 sample
% pk_rec = max(pk_rec, 1);
% pk_rec = unique(pk_rec, 'stable');
% pk_rec = sort(pk_rec);
% 
% obj.HrPeaksAfterKalman = pk_rec;
% 
% end
%

% second implement
function KalmanFilterBeats(obj)
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

    P = 10;           % Initial uncertainty
    Q = 0.5;          % Process noise (Standard deviation of beat-to-beat change)
    R_base = 5.0;     % Measurement noise (Trust in the radar peak location)

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

%MORE ADAPTIVE KALMAN FILTER:
% function KalmanFilterBeats(obj)
%     % KalmanFilterBeats (Adaptive)
%     % 1. Adaptive R: Increases when local variance of measurement is high.
%     % 2. Adaptive Q: Increases when prediction error is consistently high (Trend detection).
%     % 3. Robust Initialization & Phase Locking.
% 
%     % ---- 1. Guards & Data Prep ----
%     if isempty(obj.HrPeaks) || length(obj.HrPeaks) < 3
%         obj.HrPeaksAfterKalman = obj.HrPeaks;
%         obj.HrEstAfterKalman = [];
%         return;
%     end
% 
%     fs = obj.fs_new;
%     if isempty(fs) || fs == 0, fs = 100; end
% 
%     % Get raw data
%     raw_pk_times = unique(sort(obj.HrPeaks(:)), 'stable');
% 
%     % Calculate instantaneous HR (Measurement z)
%     ibi_meas = diff(raw_pk_times);
%     hr_meas  = 60 ./ ibi_meas; 
% 
%     % ---- 2. Initialization ----
%     % State x = [HR]
%     % Find clean start
%     valid_mask = (hr_meas > 40 & hr_meas < 200 & isfinite(hr_meas));
%     if ~any(valid_mask)
%         x = 75; 
%     else
%         x = median(hr_meas(valid_mask)); % Robust start
%     end
% 
%     P = 10;                % Initial covariance
% 
%     % Base parameters (tuning knobs)
%     Q_base = 0.5;          % Minimum process noise (steady state)
%     R_base = 5.0;          % Minimum measurement noise
% 
%     % For adaptivity: Window size for variance calculation
%     win_size = 5; 
%     recent_innovations = zeros(win_size, 1); 
% 
%     hr_hat = nan(size(hr_meas));
% 
%     % ---- 3. Adaptive Filtering Loop ----
%     for k = 1:length(hr_meas)
%         z = hr_meas(k);
% 
%         % --- A. Adaptivity Logic ---
% 
%         % 1. Adaptive R (Measurement Noise)
%         % Look at variance of recent 5 beats. 
%         % If variance is high -> signal is noisy -> Increase R (Trust Model).
%         start_idx = max(1, k - win_size);
%         local_segment = hr_meas(start_idx:k);
% 
%         % Calculate local variance (handle single point case)
%         if length(local_segment) > 1
%             local_var = var(local_segment);
%         else
%             local_var = 0;
%         end
% 
%         % Scale R: Base + factor * local_variance
%         R_adaptive = R_base + (0.5 * local_var); 
%         % Cap R to prevent filter from freezing completely
%         R_adaptive = min(R_adaptive, 100); 
% 
%         % 2. Adaptive Q (Process Noise)
%         % Calculate Innovation (Error between measurement and prediction)
%         innovation = z - x;
% 
%         % Store innovation history
%         recent_innovations = [recent_innovations(2:end); innovation];
% 
%         % If innovations are consistently large (mean error is high), 
%         % the signal is moving away (trend). Increase Q to track faster.
%         avg_innovation = mean(abs(recent_innovations));
% 
%         % Scale Q: Base + factor * average_error
%         Q_adaptive = Q_base + (0.2 * avg_innovation^2);
%         % Cap Q to prevent instability
%         Q_adaptive = min(Q_adaptive, 50);
% 
%         % --- B. Standard Kalman Steps ---
% 
%         % Predict
%         x_pred = x;
%         P_pred = P + Q_adaptive; % Use Dynamic Q
% 
%         % Gate & Update
%         is_physio = (z >= 40 && z <= 200 && isfinite(z));
% 
%         S = P_pred + R_adaptive; % Use Dynamic R
% 
%         % Gate: Mahalanobis distance
%         % If Q is high (tracking mode), we open the gate wider.
%         gate_threshold = 9 + (Q_adaptive); 
% 
%         if (innovation^2) < (gate_threshold * S) && is_physio
%              % Valid update
%              K = P_pred / S;
%              x = x_pred + K * innovation;
%              P = (1 - K) * P_pred;
%         else
%              % Outlier rejection (Keep prediction)
%              x = x_pred;
%              P = P_pred;
%         end
% 
%         hr_hat(k) = x;
%     end
% 
%     obj.HrEstAfterKalman = hr_hat;
% 
%     % ---- 4. Phase-Locked Reconstruction (Same as before) ----
%     rec_peaks = zeros(length(hr_hat)+1, 1);
%     rec_peaks(1) = raw_pk_times(1);
%     curr_time = raw_pk_times(1);
% 
%     for k = 1:length(hr_hat)
%         pred_ibi = 60 / hr_hat(k);
%         next_pred = curr_time + pred_ibi;
% 
%         % Window: +/- 20% of the beat interval
%         window = 0.2 * pred_ibi;
%         [min_err, idx] = min(abs(raw_pk_times - next_pred));
% 
%         if min_err < window
%             next_time = raw_pk_times(idx); % Snap to real peak
%         else
%             next_time = next_pred;         % Use filtered prediction
%         end
% 
%         rec_peaks(k+1) = next_time;
%         curr_time = next_time;
%     end
%     obj.HrPeaksAfterKalman = rec_peaks;
% end


%% PLOTS        
       % plotting functions       
        %plot the Hr peaks and the signals
       % function [] = plotHrpeaks(obj) % plot only the HRpeaks
       % function [] = plotBA(obj) %plot BlandAltman plot
       % function [] = plotRR(obj) % plots resipration rate plot 
       % function [] plotDashBoard(obj) 
            % plot the following 
            % 1. radar heart signal after filter(with peaks), 
            % 2. ecg refernce signal with peaks 
            % 3. Heart rate comparison between the GT and the calculated IBI. (all three are linked axis.
            % 4. bland altman plot 
        
        % ---------------------------------------------------------
        % Plotting Functions
        % ---------------------------------------------------------

        %% Plot only the HR peaks and the signals (Radar vs ECG)
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
            title (['Ground Truth-Radar Pointwise Correlation for', name]);
            lengthmax= min(length(HrToCompare),length(obj.HrGtEstAfterMedian));
            plot(obj.HrGtEstAfterMedian(1:lengthmax),HrToCompare(1:lengthmax),"DisplayName",' Median Ground truth to RADAR heart rate'...
                ,'Marker','*','LineStyle','none')
            xlabel('ground truth HR post Median filter')
            ylabel(['RADAR estimated HR of', name])
            
            plot((30:1:150),(30:1:150));
            legend('Measurements', 'Corr Axis','Location','best'); grid on; hold off;
        end       
        function h = plotHrpeaks(obj)
            h = figure('Name', 'Radar_ECG_Peaks', 'Color', 'w');
            
            % Subplot 1: Radar Signal
 
            ax(1) = subplot(2,1,1);
            hold on;
            title(sprintf('Radar HR Filters & Peak Finder - ID: %s, Scenario: %s', string(obj.ID), obj.sceneario));
            
            plot(obj.vTimeNew, obj.HrSignal*1e4, 'b', 'DisplayName', 'Filtered Radar');
            
            % Interpolate to find amplitude at peak times
            if ~isempty(obj.HrPeaks)
                peakAmps = interp1(obj.vTimeNew, obj.HrSignal, obj.HrPeaks); %changed from HrPeaks
                %peakAmps1 = interp1(obj.vTimeNew, obj.HrSignal, obj.HrPeaksAfterKalman); 
                plot(obj.HrPeaks, peakAmps*1e4, 'r*', 'MarkerSize', 8, 'DisplayName', 'Radar Peaks'); %changed from HrPeaks
                %plot(obj.HrPeaksAfterKalman, peakAmps1*1e4, 'g*', 'MarkerSize', 8, 'DisplayName', 'Radar Peaks after spike fix'); %changed from HrPeaks
            end
            
            ylabel('Amp (scaled)');
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
            
            xlabel('Time(s)'); ylabel('Amp');
            legend('show', 'Location', 'best'); grid on; hold off;
            
            linkaxes(ax, 'x');
        end

        %% Plot Bland-Altman plot
        function h = plotBA(obj, name,HrToCompare)
            h = []; % Default empty if no data
            % Check if estimates exist
            if isempty(obj.HrEst) || isempty(obj.HrGtEst)
                warning('HR Estimates are empty. Run FindRates() first.');
                return;
            end

            h = figure('Name', 'Bland_Altman_Analysis', 'Color', 'w');

            if ~isempty(obj.HrGtEst) % && ~isempty(obj.HrEstAfterKalman
                if exist('BlandAltman', 'file')
                   BlandAltman(HrToCompare,obj.HrGtEstAfterMedian, 2, 0);
                % else
                %     plot(vec_ecg, vec_radar, 'o', 'DisplayName', 'Data Points'); 
                %     xlabel('ECG'); ylabel('Radar');
                %     legend('show', 'Location', 'best');
                end
                title(sprintf('Bland-Altman of %s - ID: %s, Scenario: %s',name,string(obj.ID), obj.sceneario));
            end
                title(sprintf('Bland-Altman - ID: %s, Scenario: %s', string(obj.ID), obj.sceneario));
        end

       %% Plots respiration rate signal (Split into 2 Figures)
       %% Plot Respiration Rates (Trends)
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
        %% Plot Respiration Signals (Time Domain)
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
            xlabel('Time(s)');
            legend('show', 'Location', 'best'); grid on; hold off;          
            
            linkaxes(ax, 'x');
        end    

       %% Plot the full dashboard (Summary)
       function h = plotDashBoard(obj,name,HrToCompare) 
            h = figure('Name', 'Summary_Dashboard', 'Units', 'normalized', 'Position', [0.1 0.1 0.6 0.8], 'Color', 'w');
            
            ax_link = [];

            % 1. Radar Heart signal WITH PEAKS
            ax_link(1) = subplot(4,1,1);
            hold on;
            plot(obj.vTimeNew, obj.HrSignal, 'b', 'DisplayName', 'Radar Signal');
            if ~isempty(obj.HrPeaks)
                peakAmps = interp1(obj.vTimeNew, obj.HrSignal, obj.HrPeaks); %changed form HrPeaks
                plot(obj.HrPeaks, peakAmps, 'r*', 'MarkerSize', 8, 'DisplayName', 'Radar Peaks');
                % peakAmps = interp1(obj.vTimeNew, obj.HrSignal, obj.HrPeaksFinal); %changed form HrPeaks
                % plot(obj.HrPeaksFinal, peakAmps, 'r*', 'MarkerSize', 8, 'DisplayName', 'Radar Peaks');
            end
            title(sprintf('Radar Signal - ID: %s, Scenario: %s', string(obj.ID), obj.sceneario));
            ylabel('Amp'); legend('show', 'Location', 'best'); grid on; axis tight; hold off;

            % 2. ECG reference WITH PEAKS
            ax_link(2) = subplot(4,1,2);
            hold on;
            plot(obj.vTimeOriginal, obj.ecg_gt, 'k', 'DisplayName', 'ECG Signal');
            if ~isempty(obj.ecgPeaks)
                peakAmpsEcg = interp1(obj.vTimeOriginal, obj.ecg_gt, obj.ecgPeaks);
                plot(obj.ecgPeaks, peakAmpsEcg, 'r*', 'MarkerSize', 8, 'DisplayName', 'QRS Peaks');
            end
            title(sprintf('ECG Reference - ID: %s, Scenario: %s', string(obj.ID), obj.sceneario));
            ylabel('Amp'); legend('show', 'Location', 'best'); grid on; axis tight; hold off;

            % 3. Heart Rate Comparison (BPM vs Time)
            ax_link(3) = subplot(4,1,3);
            hold on;
            if ~isempty(obj.HrGtEst)
                time_gt_bpm = obj.ecgPeaks(2:end); 
                plot(time_gt_bpm, obj.HrGtEst, 'r.-', 'LineWidth', 1.5, 'DisplayName', 'ECG GT');
            end
            % if ~isempty(obj.HrEstFinal)
            %     time_radar_bpm = obj.HrPeaksFinal(2:end);
            %     plot(time_radar_bpm, obj.HrEstFinal, 'b.--', 'LineWidth', 1.2, 'DisplayName', 'HrEstSpikes');
            % end
              time_radar_bpm = obj.HrPeaks(2:end); %% TODO: Verify time vector! ! ! ! !!
             plot(time_radar_bpm, obj.HrEst, 'b.--', 'LineWidth', 1.2, 'DisplayName', 'HrEstSpikes');
            title(sprintf('Heart Rate (BPM) - ID: %s, Scenario: %s', string(obj.ID), obj.sceneario));
            xlabel('Time (s)'); ylabel('BPM'); legend('show', 'Location', 'best'); grid on; axis tight; hold off;

            % 4. Bland-Altman (NOT LINKED)
            subplot(4,1,4);
            if ~isempty(obj.HrGtEst) && ~isempty(obj.HrEst)
                if exist('BlandAltman', 'file')
                   BlandAltman(HrToCompare,obj.HrGtEstAfterMedian, 2, 0);
                % else
                %     plot(vec_ecg, vec_radar, 'o', 'DisplayName', 'Data Points'); 
                %     xlabel('ECG'); ylabel('Radar');
                %     legend('show', 'Location', 'best');
                end
                title(sprintf('Bland-Altman - ID: %s, Scenario: %s', string(obj.ID), obj.sceneario));
            end

            linkaxes(ax_link, 'x');            
        end
        function h = plot_examples(obj)
          % TODO - REVAMP DRAWINGS AND TITLES TO FIT THE WANTED SIGNAL. 
         h = gobjects(0); % Initialize graphics array

         % --- FIGURE 1: Peaks Comparison on HrSignal & ECG ---
         h(end+1) = figure('Name', 'Peaks_Comparison', 'Color', 'w');
         
         ax_link = [];

         % Subplot 1: Radar Signal with Initial vs Final Peaks
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
         legend('show', 'Location', 'best'); grid on; hold off;

         % Subplot 2: ECG Signal with ECG Peaks vs Radar Peaks (Matched)
         ax_link(2) = subplot(2,1,2);
         hold on;
         title('ECG Signal: GT Peaks vs Correlated Radar Peaks');
         plot(obj.vTimeOriginal, obj.ecg_gt, 'k', 'DisplayName', 'ECG Signal');
         
         if ~isempty(obj.ecgPeaks)
             amps = interp1(obj.vTimeOriginal, obj.ecg_gt, obj.ecgPeaks);
             plot(obj.ecgPeaks, amps, 'bo', 'MarkerSize', 6, 'DisplayName', 'ECG Peaks (GT)');
         end       
         % Subplot 3: Radar Signal with Initial Peaks vs GT Peaks (Matched)
         % ax_link(3) = subplot(3,1,3);
         % hold on;
         % title('Radar Signal: Radar Peaks vs Correlated GT Peaks');
         % plot(obj.vTimeNew, obj.HrSignal, 'b', 'DisplayName', 'Radar Signal');
         % 
         % if ~isempty(obj.HrPeaks)
         %    amps = interp1(obj.vTimeNew, obj.HrSignal, obj.HrPeaks);
         %    plot(obj.HrPeaks, amps, 'g*', 'MarkerSize', 6, 'DisplayName', 'Radar Peaks');
         % end
         legend('show', 'Location', 'best'); grid on; hold off;
         linkaxes(ax_link, 'x');

         % --- FIGURE 2: Heart Rate Comparisons ---
         h(end+1) = figure('Name', 'Heart_Rate_Comparison', 'Color', 'w');
         
         % Create Hr_corr_fix
         % Logic: Calculate instantaneous HR from consecutive valid matches in correlated_HrPeaks
         % Filter for 40-140 BPM range
         gt_peaks = obj.ecgPeaks;
         radar_peaks = obj.HrPeaks;
         
         is_match = radar_peaks > 0; % Valid radar peaks (not -1)
         
         % We need consecutive valid matches to compute a diff
         valid_consecutive = is_match(1:end-1) & is_match(2:end);
         
         diffs = radar_peaks(2:end) - radar_peaks(1:end-1);
         hr_vals = 60 ./ diffs; % Instantaneous HR
         
         % Create vector aligned with GT intervals (length N-1)
         Hr_corr_fix = nan(length(gt_peaks)-1, 1);
         Hr_corr_fix(valid_consecutive) = hr_vals(valid_consecutive);
         
         % Mask 40-140
         mask_range = (Hr_corr_fix >= 40) & (Hr_corr_fix <= 140);
         Hr_corr_fix(~mask_range) = NaN; 

         ax_hr = [];

         % Subplot 1
         ax_hr(1) = subplot(2,1,1); hold on; grid on;
         title(sprintf('HR Comparison: HR after filtering vs GT vs GT after median - ID: %s, Scenario: %s', string(obj.ID), obj.sceneario));
         plot(obj.HrGtEst, 'r', 'DisplayName', 'GT Est');
         plot(obj.HrEst, 'b', 'DisplayName', 'Radar Est');
         plot(obj.HrGtEstAfterMedian, 'g--','DisplayName', 'GT After Median');
         legend('show', 'Location', 'best');

         % Subplot 2
         ax_hr(2) = subplot(2,1,2); hold on; grid on;
         title('HR after filtering and preproccessing vs GT vs GT after median');
         plot(obj.HrGtEst, 'r', 'DisplayName', 'GT Est');
          if ~isempty(obj.HrEstAfterKalman)
             plot(obj.HrEstAfterKalman, 'b','DisplayName', 'Radar Est after kalman');
          end
         plot(obj.HrGtEstAfterMedian, 'g--',  'DisplayName', 'GT After Median');
         legend('show', 'Location', 'best');

         % Subplot 3
         % ax_hr(3) = subplot(4,1,3); hold on; grid on;
         % title('HR after median vs GT vs GT after median');
         % plot(obj.HrGtEst, 'r', 'DisplayName', 'GT Est');
         % if ~isempty(obj.HrEstAfterMedian)
         %     plot(obj.HrEstAfterMedian, 'b','DisplayName', 'Radar Est After Median');
         % end
         % plot(obj.HrGtEstAfterMedian, 'g','DisplayName', 'GT After Median');
         % legend('show', 'Location', 'best');
         % 
         % % Subplot 4
         % ax_hr(4) = subplot(4,1,4); hold on; grid on;
         % title(' HR (Correlated Fix) vs GT vs GT after median');
         % plot(obj.HrGtEst, 'r', 'LineWidth', 1, 'DisplayName', 'GT Est');
         % plot(obj.HrGtEstAfterMedian, 'g', 'LineWidth', 1, 'DisplayName', 'GT After Median');
         % plot(Hr_corr_fix, 'b.-', 'LineWidth', 1, 'DisplayName', 'Hr Corr Fix');
         % legend('show', 'Location', 'best');
         % xlabel('Beat Index'); ylabel('BPM');

         linkaxes(ax_hr, 'x');
     end
%% Plot Error Analysis (Corrected Titles)
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
            ylabel('Error (BPM)');
            legend('Raw Error', 'Zero Line', 'Location', 'best');
            
            % --- Subplot 2: Absolute Error vs HR (The "MAE" request) ---
            subplot(3,1,2);
            hold on; grid on;
            if ~isempty(obj.mae2Hr)
                plot(obj.mae2Hr(:,1), obj.mae2Hr(:,2), 'ko', 'MarkerFaceColor', 'g', 'MarkerSize', 4);
            end
            % Updated Title with ID and Scenario
            title(sprintf('Absolute Error (Magnitude) vs ECG HR - ID: %s, Scenario: %s', string(obj.ID), obj.sceneario));
            ylabel('|Error| (BPM)');
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
            fprintf('# Missed beats: %d\n', obj.correlation_misses);
            fprintf('# Excess beats: %d\n', obj.correlation_excess);
            fprintf('------------------------------------------------\n');
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
                
                pngFileName = fullfile(saveDir, [baseFileName, '.png']);
                figFileName = fullfile(saveDir, [baseFileName, '.fig']);

                try
                    saveas(hFig, pngFileName, 'png'); 
                    saveas(hFig, figFileName, 'fig');
                    fprintf('  -> Saved: %s\n', baseFileName);
                catch ME
                    fprintf(2, 'Warning: Could not save %s. Error: %s\n', baseFileName, ME.message);
                end
             end
        end

        %% Plot All and Optionally Save
      %% Plot All with Selection Flags

      function [] = PlotAll(obj, bsave, saveDir,name,HrToCompare, options) 
            arguments
                obj
                bsave (1,1) logical = false
                saveDir (1,1) string = 'SavedAnalysisFigures'
                name {mustBeText} = ' ' %NEW, REORDER CALLS TO PLOTALL
                HrToCompare = obj.HrEst ;
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
                h5 = obj.plotDashBoard(name,HrToCompare);
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