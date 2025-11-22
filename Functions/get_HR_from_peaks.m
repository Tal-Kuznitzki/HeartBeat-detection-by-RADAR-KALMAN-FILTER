function [vOutput] = get_HR_from_peaks(vSignal, fs)

[pks,locs,widths,proms] = findpeaks(vSignal, "MinPeakHeight",max(vSignal)*0.05);

vOutput = (fs / median(diff(locs))) * 60; % Convert to beats per minute
