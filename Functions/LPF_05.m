function [pf_05] = LPF_05(fs)

pf_05 = designfilt('lowpassfir','FilterOrder',20,...
        'PassbandFrequency',0.5,'PassbandRipple',0.05, ...
        'SampleRate',fs);


%vOutput=filtfilt(pf_05,vSignal);
