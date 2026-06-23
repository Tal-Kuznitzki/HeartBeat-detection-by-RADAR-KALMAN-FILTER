function [firL,firH,lpf5] = HRfir(fs)

% Design for a FIR HPF with a StopbandFreq 0.4 and PassbandFreq 0.6
% fs for sampling freq, and N<100
N = 100; % Filter order

StopbandFrequencyH = 0.5;
StopbandFrequencyL = 4.4;

PassbandFrequencyH=0.8;
PassbandFrequencyL=2.4;

% we pass between 0.7-3, and 0.4<Signal<4


StopbandAttenuationH=80;
StopbandAttenuationL=80;

firH = designfilt('highpassfir','StopbandFrequency',StopbandFrequencyH,...
    'PassbandFrequency',PassbandFrequencyH,'StopbandAttenuation',StopbandAttenuationH, ...
        'SampleRate',fs,'passbandRipple',0.05);

firL = designfilt('lowpassfir','PassbandFrequency',PassbandFrequencyL,...
    'StopbandFrequency',StopbandFrequencyL,'StopbandAttenuation',StopbandAttenuationL, ...
        'SampleRate',fs,'passbandRipple',0.05);

lpf5 = designfilt('lowpassfir','PassbandFrequency',5,...
    'StopbandFrequency',10,'StopbandAttenuation',StopbandAttenuationL, ...
        'SampleRate',fs,'passbandRipple',0.05);


