%% trying an improved algorithm: sort peaks and then choose the first ones

function [HR] = findHRCorrFirstPeak(vSig, fs)
%1 plot xcorr
[sigCorr,lags] = xcorr(vSig,'coeff');
positive=lags>0;
sigCorr=sigCorr(positive);
lags=lags(positive);

thresh=max(sigCorr)*0.1;
%2 find p
%plot(sigCorr);
[pks,locs,w,p] = findpeaks(sigCorr,"MinPeakDistance",fs*0.3,"MinPeakHeight",thresh);
%[sortedPks, indx] = sort(pks, "descend"); 
%sortedLoc=locs([indx]);
HR2 = lags(locs(1))/fs;

%3 calculate hr from peaks

HR = 60/HR2;


