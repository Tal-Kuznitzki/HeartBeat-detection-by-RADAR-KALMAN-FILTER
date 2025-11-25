function [HR] = findHRCorr(vSig, fs)
%1 plot xcorr
sigCorr= xcorr(vSig,vSig);
thresh=max(sigCorr)*0.01;
%2 find p
%plot(sigCorr);
[pks,locs,w,p] = findpeaks(sigCorr,"MinPeakDistance",fs*0.3,"MinPeakHeight",thresh);

%3 calculate hr from peaks
timeLags= diff(locs);
timeLags = timeLags(timeLags>fs*0.3);
timeLags = timeLags(timeLags<fs*3);
HR = (fs / median(timeLags)) * 60; % Convert to beats per minute
%HR = (fs / median(diff(locs))) * 60; % Convert to beats per minute


