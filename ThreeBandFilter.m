% this function divides a signal into below heart rate, heart rate,
% and above heart rate bands
%input: a signal, and a filter kernel size (2, 4, 8...)
%
%output: three signals, filtered and decimated. and the new sample rate of
% each signal
%
function [heartRateSignal, FsHeart, lowBandSignal, FsLow, highBandSignal, FsHigh] = ThreeBandFilter(x, Fs)

%   x  - input signal (column vector)
%   Fs - original sampling rate (Hz)
%
% Outputs:
%   heartRateSignal - band-passed 0.5–3 Hz signal (decimated)
%   FsHeart         - sampling rate of heartRateSignal (~125 Hz)
%
%   lowBandSignal   - low-frequency (<0.5 Hz) trend signal (decimated)
%   FsLow           - sampling rate of lowBandSignal (~16 Hz)
%
%   highBandSignal  - mid-range (3–50 Hz) band-passed signal (decimated)
%   FsHigh          - sampling rate of highBandSignal (~250 Hz)
%
% Filters are FIR (order ≤100) applied with zero-phase (filtfilt).
% Decimation uses anti-alias protected multistage decimation.

    arguments
        x (:,1) double
        Fs (1,1) double {mustBePositive}
    end

    x = detrend(x, 'linear');     % remove slow drift / DC bias

    % Desired effective sample rates after decimation
    FsHeartTarget = 125;    % heartbeat band target rate
    FsLowTarget   = 16;     % slow trend rate
    FsHighTarget  = 250;    % mid/high band target rate

    % Compute integer decimation factors
    Rheart = max(1, round(Fs / FsHeartTarget));
    FsHeart = Fs / Rheart;

    Rlow = max(1, round(Fs / FsLowTarget));
    FsLow = Fs / Rlow;

    Rhigh = max(1, round(Fs / FsHighTarget));
    FsHigh = Fs / Rhigh;

    % Apply multistage decimation with proper anti-alias filtering
    heartBase = decimate_multistage(x, Rheart);   % -> FsHeart
    lowBase   = decimate_multistage(x, Rlow);     % -> FsLow
    highBase  = decimate_multistage(x, Rhigh);    % -> FsHigh

    % Nyquist limits after decimation
    nyHeart = FsHeart/2;
    nyLow   = FsLow/2;
    nyHigh  = FsHigh/2;

    % ----------- (1) Heart-rate band: 0.5–3 Hz -----------
    pass1 = [0.5 3];
    stop1 = [0.3 4.0];  % widened for transition

    dHeart = designfilt('bandpassfir', ...
        'FilterOrder', 100, ...
        'StopbandFrequency1', stop1(1), 'PassbandFrequency1', pass1(1), ...
        'PassbandFrequency2', pass1(2), 'StopbandFrequency2', stop1(2), ...
        'StopbandAttenuation1', 60, 'StopbandAttenuation2', 60, ...
        'PassbandRipple', 0.5, ...
        'SampleRate', FsHeart);

    heartRateSignal = filtfilt(dHeart, heartBase);

    % ----------- (2) Low-frequency band: < 0.5 Hz -----------
    pass2 = min(0.45, nyLow * 0.45);
    stop2 = min(0.70, nyLow * 0.9);

    dLow = designfilt('lowpassfir', ...
        'FilterOrder', 100, ...
        'PassbandFrequency', pass2, ...
        'StopbandFrequency', stop2, ...
        'PassbandRipple', 0.5, ...
        'StopbandAttenuation', 70, ...
        'SampleRate', FsLow);

    lowBandSignal = filtfilt(dLow, lowBase);

    % ----------- (3) Mid-high band: 3–50 Hz -----------
    pass3 = [3 50];
    stop3 = [2 60];

    if pass3(2) >= nyHigh
        error('After decimation FsHigh=%.1f Hz => Nyquist=%.1f < 50 Hz. Choose different decimation.', FsHigh, nyHigh);
    end

    dHigh = designfilt('bandpassfir', ...
        'FilterOrder', 100, ...
        'StopbandFrequency1', stop3(1), 'PassbandFrequency1', pass3(1), ...
        'PassbandFrequency2', pass3(2), 'StopbandFrequency2', stop3(2), ...
        'StopbandAttenuation1', 60, 'StopbandAttenuation2', 70, ...
        'PassbandRipple', 0.5, ...
        'SampleRate', FsHigh);

    highBandSignal = filtfilt(dHigh, highBase);
end


% ---------- Helper: decimate in small safe stages ----------
function y = decimate_multistage(x, R)
    if R < 2
        y = x;
        return;
    end
    primesSmall = [2 3 5];
    f = factor(R);
    leftover = prod(f(~ismember(f,primesSmall)));
    f = f(ismember(f,primesSmall));
    if isempty(f), f = 1; end
    stages = [f(:).' leftover];
    stages(stages == 1) = [];

    y = x;
    for k = 1:numel(stages)
        y = decimate(y, stages(k), 'fir');
    end
end

