radar_fs=fs_radar;
vRadarPhase = atan2(radar_q,radar_i);
vRadarAmp=tfm_ecg; %sqrt(radar_q.^2 + radar_i.^2);
vRadarAmp = vRadarAmp;
vRadar_dist = radar_dist;

radarFFt1=fftshift(fft(vRadarAmp));
maxN=length(vRadar_dist);
N2=maxN/2;
radar2_padded=[vRadar_dist(1:N2); zeros(N2,1)];
radarFFt2=fftshift(fft(phase_compensated));
% moving average on radarAmp:
vNormRadarAmped=radar_dist;
for i=1000:1000:maxN
    vNormRadarAmped(i-999:i)=vNormRadarAmped(i-999:i)-mean(vNormRadarAmped(i-999:i));
end
radarFFt3=fftshift(fft(vNormRadarAmped));

vNormXCorr=xcorr(vNormRadarAmped,vNormRadarAmped);


xFreq1=(-0.5:1/maxN:0.5-1/maxN)*radar_fs;

xFreq2=(-0.5:1/N2:0.5-1/N2)*radar_fs;

figure(1) %fft of full and half signal
subplot(4,1,1)
plot(xFreq1,abs(radarFFt3));

subplot(4,1,2)
plot(vNormRadarAmped);

subplot(4,1,3)
plot(xLags,abs(vNormXCorr));
title('corr for normalized radar_dist')

subplot(4,1,4)
plot(vRadarPhase);

hold off;

vAmpCorr= xcorr(vNormRadarAmped,vNormRadarAmped);
xLags=(-1+1/maxN:1/(maxN):1-1/maxN)*(maxN/radar_fs);
figure(2);
plot(xLags,vAmpCorr);


RadarFs100=decimate(phase_compensated,20);

filtRadarFs100=conv(RadarFs100,LPF3Hz100Fs,'same');

fftPostLPF100=fftshift(fft(filtRadarFs100));
xFreq3=(-0.5:(20/maxN):0.5-(20/maxN))*radar_fs/20;

%now HPF to find BPM
filteredFs20=decimate(filtRadarFs100,5);
filt2HPF20=conv(filteredFs20,HPF05Hz20Fs,'same');

xFreq4=(-0.5:(100/maxN):0.5-(100/maxN))*radar_fs/100;

figure(3)
title('decimate 20, lpf 3Hz')
subplot(2,1,1)
plot(filtRadarFs100);
subplot(2,1,2)
plot(xFreq3,abs(fftPostLPF100));

filt2FFT=fftshift(fft(filt2HPF20));

figure(4)
title('decimate 100, lpf and hpf')
subplot(2,1,1)
plot(filt2HPF20);
subplot(2,1,2)
plot(xFreq4,abs(filt2FFT));

%stft on phase
win = hann(8192, "periodic");
noverlap = 4096;   % 50% overlap (typical)
phase_decimated=decimate(phase_compensated,10);
phaseSTFT1=stft(phase_decimated, radar_fs, "Window", win);

figure(5)
mesh(abs(phaseSTFT1));

