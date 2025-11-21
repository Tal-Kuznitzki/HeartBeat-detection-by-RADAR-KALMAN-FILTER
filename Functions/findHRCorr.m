function [HR] = findHRCorr(vSig, fs)
%1 plot xcorr
sigCorr= xcorr(vSig,vSig);

%2 find p
%plot(sigCorr);
[pks,locs,w,p] = findpeaks(sigCorr,"MinPeakDistance",fs*0.6);

%3 calculate hr from peaks
HR = (fs / mean(diff(locs))) * 60; % Convert to beats per minute


