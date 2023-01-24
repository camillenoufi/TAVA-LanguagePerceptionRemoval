%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Generation of TAVA 'Auditory' Stimuli from EGG + VT   %
%            via MATLAB Implementation                 %
%                                                      %
% Author: Camille Noufi                      4/26/2021 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% interface function, run from matlab command line
function [input_path, output_path]  = generate_stimulus(input_path)
[path,name,ext] = fileparts(input_path);
output_path = fullfile(path, [name,'-music',ext]);

% Import external scripts & helper functions.  
% change this if directory structure changes
addpath ('./helpers')
addpath('./../voicebox');
addpath('./../GFM-IAIF-master');

% FLAGS TO GENERATE BASELINE COMPARISONS:
lpf_flag = 1;
lpc_source = 1;

% define STFT processing params
fs_lpc = 16e3; % fs<=16khz needed for accurate LPC
% AM params
wlen = 512;
hop = wlen/8;

% initialize LPC / filter params
lpc_params.p = 2 + fs_lpc/2/1e3; %LPC order
lpc_params.p_glot = 3;

%% Main Wrapper:
% Perform analysis of speech and synthesized of paired "music" stimuli 

% ~~~~~~~~~~~~~~~~~ get EGG & speech mono signal ~~~~~~~~~~~~~~~~~~~~~~~~~
[signal, fs, y_name] = LoadSpeechSignal(input_path);
g = NormalizeLoudness(signal(:,1),fs); %EGG
s = NormalizeLoudness(signal(:,2),fs); %Speech

% ~~~~~~~~~~~~~~~~~~~~ VT Shaping ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% loop through voiced frames to find average vocal tract transfer function
g_vt = ApplyAverageLPCTransferFunction(s, g, fs_lpc, fs, lpc_params);

% ~~~~~~~~~~~~~~~~~~~ Amplitude Modulation ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% fine cross-filter to apply amplitude modulation of EGG by speech envelope,
% same windowing as filtering

[g_am,~,~] = ModulateAmplitude(g_vt,s,fs,wlen,hop);
g_am(isnan(g_am)) = 0; %set NaNs to 0;
g_am = NormalizeLoudness(g_am',fs);

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ save ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
sigs = {g(:), s(:), g_am(:)}; 
sigs = EqualizeLengths(sigs); %trim to be same length for multichannel comparison
audiowrite(output_path, sigs{3}, fs);
end

%% Helper Functions
   
%load speech signal
function [y,fs,name] = LoadSpeechSignal(input_filepath)
    % read file
    [~,name,~] = fileparts(input_filepath);
    [y, fs] = audioread(input_filepath);
    disp(name);
end

% determine window gain for short-time operations to ensure COLA
function winGain = GetWindowGain(win,N,M,R)
    % get window-scaling factor for COLA
    winGain = zeros(N,1); 
    for s=0:R:3*M
        ndx = s+1:s+M; %current window location
        winGain(ndx) = winGain(ndx) + win;  %overlap-add energy window
    end
    winGain = max(winGain);
end

function y = NormalizeLoudness(x,fs)
    [loudness, ~] = integratedLoudness(x,fs);
    target = -23; %EBU R 128 standard
    gaindB = target - loudness;
    gain = 10^(gaindB/20);
    y = x.*gain;
end

function y = NormalizePeak(x)
    maxval = max(abs(x));
    y = 0.99*x/maxval;
end

function g_vt = ApplyAverageLPCTransferFunction(s, g, fs, fs_original, lpc_params)

    %unpack
    p = lpc_params.p;
    p_glot = lpc_params.p_glot;

    %define window-len M and hop size R for COLA method
    wlen = 2*fs*10e-3; %20ms window
    M = wlen+1;
    R = (M-1)/2;
    win = hamming(M); win(1) = win(1)/2; win(M) = win(M)/2; %to ensure COLA

    % downsample speech signal to work better for LPC
    s = resample(s,fs,fs_original);

    % loop setup
    H = zeros(wlen,1);
    [~,~,~,vad]=v_activlev(s,fs);
    sVoiced = s(vad);
    N = numel(sVoiced);
    winGainCOLA = GetWindowGain(win,N,M,R);

    % loop to get LPC TF per frame
    n = 1; %frame idx
    for j=0:R:(N-M)
        ndx = j+1:j+M; %current window location
        sn = (winGainCOLA*win).*sVoiced(ndx); %get windowed speech sample
        [Av_n,~,~] = gfmiaif(sn,p,p_glot, 0.99, win); %get LPC coeffs for vt, gs, and lr
        [Hn,F] = freqz(1,Av_n,wlen,fs); %get freq-resp amplitudes and f-locations based on LPC coefficients
        if isnan(Hn)
            continue %pass disp(n); disp('NaN calculated') %disp(j+1); disp(j+M)
        else
            H = H + abs(Hn); %add to sum
        end
        n = n+1; %increment frame number
    end
    H = H/(n-1); %avg
    [avgIR,avgTF] = MinPhaseFIR(fs,F,H,256); %create minimum-phase filter, right now a 1024pt filter before upsampling

    %resample in the time domain to EGG's original sampling rate
    avgIR_resampled = interp(avgIR,fs_original/fs);
    avgIR_resampled = NormalizePeak(avgIR_resampled); %max filter amplitude at 1

%     % BASELINE:
%     if lpc_source
%         
%     end
        
    % %convolve filter with EGG signal at original sampling rate
    g_vt = conv(g,avgIR_resampled,'same');

end

% minimum phase FIR filter given amplitudes and frequencies
function [IR,mag_fullspectrum] = MinPhaseFIR(fs, fc, A, taps)

    fir.freq_vec = fc;
    fir.magnitude = A;

    filter_len = 2^(nextpow2(fs/2)); %filter length must be an even number

    freq_vec = fs*(0:filter_len-1)./filter_len;

    mag_resampled = interp1(fir.freq_vec, fir.magnitude, freq_vec(1:(filter_len/2)+1)); %interpolate between fc samples through nyquist
    %replace nans if they exist
    mag_resampled = smooth(fillmissing(mag_resampled,'linear','SamplePoints',freq_vec(1:length(mag_resampled))),'moving');
    nyquist_upper_half = flipud(mag_resampled);

    mag_fullspectrum = [mag_resampled(1:filter_len/2); nyquist_upper_half(1:filter_len/2)];
    mag_fullspectrum = mag_fullspectrum + 1e-6; %add eps to prevent log(0) next line

    % calculate phase through hilbert transform
    phase = imag(-hilbert(log(mag_fullspectrum)));
    complex_spectrum = mag_fullspectrum.*exp(1i*phase);
    IR = real(ifft(complex_spectrum));

    % not normalized
    IR = IR(1:taps);
end

% amplitude modulation through STFT cross-filtering of total energy across
% frequency bands
function [z,Z_stft,F] = ModulateAmplitude(x,y,fs, wlen, hop)
    
    nfft = wlen;

    % perform time-frequency analysis
    [X_stft, ~] = stft(x, wlen, hop, nfft, fs); %carrier
    [Y_stft, F] = stft(y, wlen, hop, nfft, fs); %modulator

    % extract envelope of the modulating signal y
    Y_energy = sum(abs(Y_stft),1);
    Y_energy = Y_energy/max(Y_energy);
    % extract envelope of the carrier signal x
    X_energy = sum(abs(X_stft),1);
    X_energy = X_energy/max(X_energy);
    % calculate scaling ratio between modulator and carrier
    R_vec = Y_energy./X_energy;
    % attenuate to EGG level if no EGG energy but speech energy (likely all consonants)
    R_vec(R_vec > 2) = 1;
    
    % convert vector to diagonal right-multiply matrix to affect columns of
    % X
    R = diag(R_vec);
    
    % apply ratio to all frames in X to scale to modulator Y's energy
    Z_stft = X_stft*R;
    
    %back to time domain
    z = istft(Z_stft, wlen, hop, nfft, fs);   
end

%equalize lengths of mono audio signals so they can be combined into
%multichannel output wav file
function sigs = EqualizeLengths(sigs)
   
    lens = zeros(length(sigs),1);
    for k = 1:length(sigs)
        lens(k) = min(length(sigs{k}));
    end
    minlen = min(lens);

    for k = 1:length(sigs)
        if length(sigs{k}) > minlen
            tmp = sigs{k};
            sigs{k} = tmp(1:minlen);
        end
    end

end