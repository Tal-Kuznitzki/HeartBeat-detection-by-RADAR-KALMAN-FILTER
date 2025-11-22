function [vOutput] = get_HR_from_peaks(vSignal, fs)
thresh=mean(abs(vSignal))*0.1;

[pks,locs,widths,proms] = findpeaks(vSignal, "MinPeakHeight",thresh,'MinPeakDistance',0.33*fs);

vOutput = (fs / median(diff(locs))) * 60; % Convert to beats per minute
