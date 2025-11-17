function [vOutput] = BPF_005_05(vSignal, fs)

pf_05 = designfilt('lowpassfir','FilterOrder',40,...
        'CutoffFrequency1',0.05,'CutoffFrequency2',0.5, ...
        'SampleRate',fs);

vOutput=filtfilt(pf_05,vSignal);
