function [pf_05] = LPF_05(fs)

pf_05 = designfilt('lowpassfir','FilterOrder',100,...
        'PassbandFrequency',0.5, ...
        'SampleRate',fs);

%'PassbandRipple',0.05,
%vOutput=filtfilt(pf_05,vSignal);
