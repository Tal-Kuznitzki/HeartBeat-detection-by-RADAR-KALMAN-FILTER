function [rr_bpm, t_rr, rr_interp] = RespRateFromPeaks(peakLocs, winSec, fs)
% RespRateFromPeaks
%
% Estimates respiration rate from detected respiration peaks.
%
% Method:
%   - Use a sliding window, default 10 seconds
%   - Window center moves 1 second at a time
%   - For every second, count how many peaks fall inside the window
%   - Convert peak count to breaths per minute
%   - Interpolate over the full time axis to smooth/complete the estimate
%
% Inputs:
%   peakLocs - respiration peak locations
%              default: in seconds
%              if fs is provided: in samples
%
%   winSec   - window length in seconds
%              default: 10
%
%   fs       - optional sampling frequency
%              use only if peakLocs are sample indices
%
% Outputs:
%   rr_bpm    - respiration rate estimate at 1-second resolution [BPM]
%   t_rr      - time vector for rr_bpm [sec]
%   rr_interp - interpolated/smoothed respiration rate [BPM]

    if nargin < 2 || isempty(winSec)
        winSec = 10;
    end

    if nargin < 3
        fs = [];
    end

    peakLocs = peakLocs(:);

    if isempty(peakLocs)
        rr_bpm = [];
        t_rr = [];
        rr_interp = [];
        return;
    end

    % ------------------------------------------------------------
    % Convert peaks to seconds if they were given in samples
    % ------------------------------------------------------------
    if ~isempty(fs)
        peakTimes = peakLocs ./ fs;
    else
        peakTimes = peakLocs;
    end

    peakTimes = peakTimes(isfinite(peakTimes));
    peakTimes = sort(peakTimes(:));

    if isempty(peakTimes)
        rr_bpm = [];
        t_rr = [];
        rr_interp = [];
        return;
    end

    % ------------------------------------------------------------
    % Build 1-second time grid
    % ------------------------------------------------------------
    tStart = floor(min(peakTimes));
    tEnd   = ceil(max(peakTimes));

    t_rr = (tStart:1:tEnd).';

    if numel(t_rr) < 2
        rr_bpm = NaN(size(t_rr));
        rr_interp = rr_bpm;
        return;
    end

    halfWin = winSec / 2;

    rr_bpm = nan(size(t_rr));

    % ------------------------------------------------------------
    % Sliding window peak counting
    % ------------------------------------------------------------
    for k = 1:numel(t_rr)

        tCenter = t_rr(k);

        t1 = tCenter - halfWin;
        t2 = tCenter + halfWin;

        nPeaks = sum(peakTimes >= t1 & peakTimes < t2);

        % breaths / sec -> breaths / minute
        rr_bpm(k) = 60 * nPeaks / winSec;
    end

    % ------------------------------------------------------------
    % Interpolate / smooth over the total time
    % ------------------------------------------------------------
    validMask = isfinite(rr_bpm);

    if sum(validMask) >= 2
        rr_interp = interp1( ...
            t_rr(validMask), ...
            rr_bpm(validMask), ...
            t_rr, ...
            'pchip', ...
            'extrap');
    else
        rr_interp = rr_bpm;
    end

    % Optional mild smoothing to avoid staircase behavior
    % Smooth only non-zero points, because zero means breath-hold / apnea
    smoothWin = max(3, round(winSec / 2));
    
    nonZeroMask = isfinite(rr_interp) & rr_interp > 0;
    
    rr_smooth = rr_interp;
    
    if sum(nonZeroMask) >= smoothWin
        rr_smooth(nonZeroMask) = movmean( ...
            rr_interp(nonZeroMask), ...
            smoothWin, ...
            'Endpoints', 'shrink');
    end
    
    rr_interp = rr_smooth;
    
end