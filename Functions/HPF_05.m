function [vOutput] = HPF_05(fs)

d = designfilt('highpassiir',...
        'StopbandFrequency', 0.4,...
        'PassbandFrequency',0.7,'StopbandAttenuation',60, ...
        'SampleRate',fs,'DesignMethod',       'ellip');

[b,a] = tf(d);   % returns numerator b and denominator a
b=b(:);
a=a(:);

% assume b,a are from tf(d)
%fprintf('max abs coeff b=%g, a=%g\n', max(abs(b)), max(abs(a)));
p = roots(a);
%fprintf('max pole mag = %g\n', max(abs(p)));

%%%%y = filter(b, a, vSignal);

vOutput =d; %just return the HPF
%%%%vOutput=y;
