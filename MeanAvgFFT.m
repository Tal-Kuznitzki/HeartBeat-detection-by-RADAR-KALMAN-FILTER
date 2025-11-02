function [maxFound] = MeanAvgFFT(vSignal, Fs)
% we want to see if we can deduce the BPM by taking half of
% the first harmonic in a moving average version of radar_dist

minFreq=0.5; %change moving average size. it is sort of a low stop filter
normWindow=minFreq*Fs-1;
vNorm=vSignal;
N=length(vSignal);
vNorm(1:normWindow)=vNorm(1:normWindow)-mean(vNorm(1:normWindow));
for i=normWindow+1:1:N
    vNorm(i)=vNorm(i)-mean(vNorm(i-normWindow:i));
end
% now vNorm is a vector centered on 0, relatively. we lost any frequency 
%below a certain threshold, determined by the window chosen.

vFftNorm= (fft(vNorm));
vAbsFftNorm=abs(vFftNorm);
figure(1)
plot(vNorm);
[value, maxFound]=max(vAbsFftNorm(normWindow:N));
maxFound=(maxFound/2)*1/(N); %supposedly we got double the GT, lets check it

