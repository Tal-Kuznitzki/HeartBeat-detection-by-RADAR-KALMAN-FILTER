function [vOutput] = get_HR_from_peaks(vSignal, fs)

[pks,locs,widths,proms] = findpeaks(vSignal, fs,"MinPeakHeight",0.5*1000);

vOutput=mean(diff(locs));
