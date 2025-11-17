function [vOutput] = HPF_05(vSignal, fs)

pf_05 = designfilt('highpassfir','FilterOrder',30,...
        'PassbandFrequency',0.5,'PassbandRipple',0.01, ...
        'SampleRate',fs);

vOutput=filtfilt(pf_05,vSignal);
