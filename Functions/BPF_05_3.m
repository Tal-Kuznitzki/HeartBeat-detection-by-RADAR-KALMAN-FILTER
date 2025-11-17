function [vOutput] = BPF_05_3(vSignal, fs)

pf_05 = designfilt('lowpassfir','FilterOrder',40,...
        'CutoffFrequency1',0.5,'CutoffFrequency2',3, ...
        'SampleRate',fs);

vOutput=filtfilt(pf_05,vSignal);
