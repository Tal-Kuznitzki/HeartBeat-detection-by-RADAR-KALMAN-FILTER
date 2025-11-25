% this function finds the HR with either peaks or correlation:
%input:
% signal
% sample freq
% method - use xCorr, find peaks  2 for QRS gt function
% step- step size in seconds

%output:
% a vector of hr every step
% a fitting time vector


function [hrVec, tVec] = movingWindowHR(x, fs,windowSec,step, method)
arguments
    
    x  {mustBeNumeric}    % accept row or column, normalize below
    fs (1,1) {mustBePositive, mustBeFinite}
    windowSec (1,1) {mustBePositive, mustBeFinite}
    step {mustBePositive}
    method (1,1) string {mustBeMember(method, ["corr","peaks","corrmax","fft","qrs"])} = "corr"
   
end

% Ensure column vector
x = x(:);

% Samples per window and hop (1 second update)
winSamps = round(windowSec * fs);
hopSamps = round(step * fs);

if winSamps < 1
    error('windowSec * fs must be >= 1 sample.');
end

N = numel(x);
if N < winSamps
    % If signal shorter than window, run single window on padded/truncated data
    pad = zeros(winSamps - N, 1);
    xP = [x; pad];
    hr = findHRCorr(xP, fs);
    hrVec = hr;
    tCenter = (N/2) / fs;  % approximate center time of original signal
    tVec = tCenter;
    return
end

% Start indices for each window: include final window that ends at N
startIdx = 1:hopSamps:(N - winSamps + 1);
if startIdx(end) + winSamps - 1 < N
    % add final window aligned to the end
    startIdx(end+1) = N - winSamps + 1;
end

nWins = numel(startIdx);
hrVec = NaN(nWins, 1);
tVec  = NaN(nWins, 1);

for k = 1:nWins
    
    s = startIdx(k);
    e = s + winSamps - 1;
    winSig = x(s:e);

    % call provided HR estimator
    try
        if(method=="corr")
            hr = findHRCorr(winSig, fs);
        elseif(method=="peaks")
            hr = get_HR_from_peaks(winSig, fs);
        elseif(method=="qrs")
            [~,hr,~] = pan_tompkin(winSig, fs,0);
            hr=(fs/median(diff(hr))) * 60;
        elseif(method=="fft")
            hr = findHRFft(winSig, fs);
        else %method is "corrmax"
            hr = findHRCorrFirstPeak(winSig,fs);
        end
    catch
        hr = NaN;
    end

    hrVec(k) = hr;
    % time stamp = center of the window (in seconds)
    tCenter = (s + e) / 2 / fs;
    tVec(k) = tCenter;
end
end
