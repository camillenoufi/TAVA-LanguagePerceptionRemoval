%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Generation of TAVA 'Auditory' Stimuli from EGG + VT   %
%            via MATLAB Implementation                 %
%                                                      %
% Author: Camille Noufi                      4/26/2021 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear variables;

% Set up file paths;
dataset_path = './../../DanStims/'; % change this to dataset folder
audio_path = fullfile(dataset_path,'Last4subjects'); % folder containing subject speech/EGG audio files each in their own subdirectory

% Output audio files will be stored in a new subfolder 'output' within each
% subject's individual folder

% Import external scripts & helper functions.  
% change this if directory structure changes
addpath ('./helpers')
addpath('./../voicebox');
addpath('./../GFM-IAIF-master');

% define STFT processing params
fs_lpc = 16e3; % fs<=16khz needed for accurate LPC
% AM params
wlen = 512;
hop = wlen/8;

% initialize LPC / filter params
lpc_params.p = 2 + fs_lpc/1e3; %LPC order
lpc_params.p_glot = 3;

%% Main Wrapper:
% Perform analysis of speech and synthesis of paired "music" 
% for each audio stimulus 
%       in each subject's folder

subjects_dir = dir(audio_path);
wb = waitbar(0,'');
macOS_start_idx = 4; %start at 4 to ignore ., .. and DS_store for Mac OS, change this to 1 if this is causing problems
K = length(subjects_dir)-(macOS_start_idx-1); %num subjects

% subject-level loop
for k = macOS_start_idx:length(subjects_dir)
    tic;
    waitbar((k-3)/K,wb,['Creating music stimuli per subject. Processing... subject: ',num2str(k-3),'/',num2str(K)]);
    
    % file management
    stimulus_folder = subjects_dir(k).name;
    output_folder = fullfile(audio_path,stimulus_folder,'output');
    if ~exist(output_folder, 'dir')
        mkdir(output_folder);
    end
    disp([num2str(k-(macOS_start_idx-1)), ': ', stimulus_folder]);
    
    % stimulus-level loop
    stimulus_dir = dir(fullfile(audio_path,stimulus_folder,'*.wav'));
    for i=1:length(stimulus_dir)

        % ~~~~~~~~~~~~~~~~~ get EGG & speech mono signal ~~~~~~~~~~~~~~~~~~~~~~~~~
        [signal, fs, y_name] = LoadSpeechSignal(stimulus_dir(i));
        g = NormalizeLoudness(signal(:,1),fs); %EGG
        s = NormalizeLoudness(signal(:,2),fs); %Speech

        % ~~~~~~~~~~~~~~~~~~~~ VT Shaping ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        % loop through voiced frames to find average vocal tract transfer function
        g_vt = ApplyAverageLPCTransferFunction(s, g, fs_lpc, fs, lpc_params);

        % ~~~~~~~~~~~~~~~~~~~ Amplitude Modulation ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        % fine cross-filter to apply amplitude modulation of EGG by speech envelope,
        % same windowing as filtering

        [g_am,G_stft,F] = ModulateAmplitude(g_vt,s,fs,wlen,hop);
        g_am(isnan(g_am)) = 0; %set NaNs to 0;
        g_am = NormalizeLoudness(g_am',fs);

        % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ save ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        % 1 - EGG
        % 2 - Speech
        % 3 - EGG with avg VT
        % 4 - EGG with avg VT shaping + AM
        
        sigs = {g(:), s(:), g_am(:), g_vt(:)}; 
        sigs = EqualizeLengths(sigs); %trim to be same length for multichannel audiowrite

        output_filepath = fullfile(output_folder,[y_name,'_music.wav']);
        % add channels to multichannel mix
        mix = [];
        for c = 1:(length(sigs)-1) %leave off intermediate VT-only signal
            mix = [mix, sigs{c}]; %create multichannel matrix
        end
        audiowrite(output_filepath, mix, fs);

    end
    toc;
end
close(wb);

%% Helper Functions

function ManageInput(varargin)
 if nargin > 1
   input = varargin{1};
   pathstr = cd();
 else
   [filename, pathstr] = uigetfile('Which file?');
 end
 fullname = fullfile(pathstr, filename);
end
   
%load speech signal
function [y,fs,name] = LoadSpeechSignal(entry)
    % read file
    [~,name,~] = fileparts(entry.name);
    [y, fs] = audioread(fullfile(entry.folder,entry.name));
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
    H = H + abs(Hn); %add to sum
    n = n+1; %increment frame number
end
H = H/(n-1); %avg
[avgIR,avgTF] = mpFIR(fs,F,H); %create minimum-phase filter, right now a 1024pt filter before upsampling

%resample in the time domain to EGG's original sampling rate
avgIR_resampled = interp(avgIR,fs_original/fs);
avgIR_resampled = NormalizePeak(avgIR_resampled); %max filter amplitude at 1

% %convolve filter with EGG signal at original sampling rate
g_vt = conv(g,avgIR_resampled,'same');

end

% minimum phase FIR filter given amplitudes and frequencies
function [IR,totalMag] = mpFIR(fs, fc, amps)

fir.freqVect = fc;
fir.magnitude = amps;

filtLen = 2^(nextpow2(fs/2)); %filter length must be an even number

freqVect = fs*(0:filtLen-1)./filtLen;

resampleMag = interp1(fir.freqVect, fir.magnitude, freqVect(1:(filtLen/2)+1)); %interpolate between fc samples through nyquist
%replace nans if they exist
resampleMag = smooth(fillmissing(resampleMag,'linear','SamplePoints',freqVect(1:length(resampleMag))),'moving');
resampleNyquist = flipud(resampleMag);

totalMag = [resampleMag(1:filtLen/2); resampleNyquist(1:filtLen/2)];
%totalMag = db2mag(totalMag);
totalMag = totalMag + 1e-6;

X = -hilbert(log(totalMag));
phase = imag(X);
complexSpec = totalMag.*exp(1i*phase);
IR = real(ifft(complexSpec));

%IR = IR./max(abs(IR));
IR = IR(1:256);

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