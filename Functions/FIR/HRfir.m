function [vSigOut] = HRfir(vSig,fs)

% Design a FIR HPF with a StopbandFreq 0.4 and PassbandFreq 0.6
% fs for sampling freq, and N<100
N = 100; % Filter order
firH = designfilt('highpassfir','StopbandFrequency',0.4,...
    'PassbandFrequency',0.7,'StopbandAttenuation',60, ...
        'SampleRate',fs);

firL = designfilt('lowpassfir','PassbandFrequency',3,...
    'StopbandFrequency',4,'StopbandAttenuation',80, ...
        'SampleRate',fs);

figure; title('HPF and LPF');
subplot(2,1,1);
plot(abs(fftshift(fft(firH.Numerator))));
subplot(2,1,2);
plot(abs(fftshift(fft(firL.Numerator))));


vSigOut= filtfilt(firH,vSig);
vSigOut = filtfilt(firL,vSigOut);