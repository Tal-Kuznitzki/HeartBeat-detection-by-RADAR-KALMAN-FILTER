function [firL,firH] = HRfir(fs)

% Design a FIR HPF with a StopbandFreq 0.4 and PassbandFreq 0.6
% fs for sampling freq, and N<100
N = 100; % Filter order

StopbandFrequencyH = 0.5;
StopbandFrequencyL = 3.7;

PassbandFrequencyH=0.8;
PassbandFrequencyL=1.7;

% we pass between 0.7-3, and 0.4<Signal<4


StopbandAttenuationH=80;
StopbandAttenuationL=80;

firH = designfilt('highpassfir','StopbandFrequency',StopbandFrequencyH,...
    'PassbandFrequency',PassbandFrequencyH,'StopbandAttenuation',StopbandAttenuationH, ...
        'SampleRate',fs,'passbandRipple',0.05);

firL = designfilt('lowpassfir','PassbandFrequency',PassbandFrequencyL,...
    'StopbandFrequency',StopbandFrequencyL,'StopbandAttenuation',StopbandAttenuationL, ...
        'SampleRate',fs,'passbandRipple',0.05);



% nameH = sprintf('FIR HPF freq response stopband Freq %.1f [Hz] Passband Freq %.1f [Hz] Stopband Attenuation %d dB',...
%     StopbandFrequencyH,PassbandFrequencyH,StopbandAttenuationH);
% nameL = sprintf('FIR LPF freq response stopband Freq %.1f [Hz] Passband Freq %.1f [Hz] Stopband Attenuation %d dB',...
%     StopbandFrequencyL,PassbandFrequencyL,StopbandAttenuationL);

%{
figure; 
set(gcf, 'Name', 'HRfir_Internal_Plot'); 
subplot(2,1,1);
plot(abs(fftshift(fft(firH.Numerator))), 'DisplayName','HPF');
title(nameH); 
legend('show');
grid on;

% Subplot 2: LPF (Low Pass Filter)
subplot(2,1,2);
plot(abs(fftshift(fft(firL.Numerator))), 'DisplayName','LPF');
title(nameL); 
legend('show');
xlabel('Frequency Bin (Shifted)');
grid on;
%}
% 
% vSigOut= filtfilt(firH,vSig);
% vSigOut = filtfilt(firL,vSigOut);