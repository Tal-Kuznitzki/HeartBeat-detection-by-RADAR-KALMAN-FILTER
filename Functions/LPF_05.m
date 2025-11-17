function [vOutput] = LPF_05(vSignal, fs)

pf_05 = designfilt('lowpassfir','FilterOrder',20,...
        'PassbandFrequency',0.5,'PassbandRipple',0.05, ...
        'SampleRate',fs);

vOutput=filtfilt(pf_05,vSignal);
