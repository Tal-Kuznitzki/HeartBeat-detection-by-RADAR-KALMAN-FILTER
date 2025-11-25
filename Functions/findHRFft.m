%% trying to find the HR using FFT 
function HR = findHRFft(vSig, fs)

    vSig = vSig(:);

    N = length(vSig);
    Nfft = 2^nextpow2(N);

    V = fft(vSig, Nfft);
    vFreq = (0:Nfft-1)*(fs/Nfft);

    % --- Search only in HR band ---
    idx = (vFreq >= 0.5 & vFreq <= 3);
    mag = abs(V(idx));
    fband = vFreq(idx);

    % --- Dominant frequency ---
    [~, i1] = max(mag);
    f1 = fband(i1);

    % --- Harmonic check ---
    f_half = f1/2;
    tol = 0.05;

    idx2 = find(abs(fband - f_half) < tol);

    if ~isempty(idx2)
        if mag(idx2) > 0.5 * mag(i1)
            f1 = f_half;
        end
    end

    HR = 60 * f1;
    if (isnan(HR))
        HR=0;
    end

end
