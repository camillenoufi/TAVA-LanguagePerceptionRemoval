% batch process generate_stimulus.m over directory

%file list
dir_path = '/Users/camillenoufi/Downloads/TrainDataClean/';
wav_dir = dir(fullfile(dir_path,'*.wav'));
output_dir = '/Users/camillenoufi/Library/CloudStorage/GoogleDrive-cnoufi@stanford.edu/My Drive/Research/TAVA/TAVA-NN/data/audio';
if ~exist(output_dir, 'dir')
    mkdir(output_dir)
end

%iterate over all files in wav_dir
f = waitbar(0, 'Starting');
n = length(wav_dir);
for i= 1:n
    % load and process speech+EGG files --> tEGG files
    file_path = fullfile(dir_path,wav_dir(i).name);
    generate_stimulus(file_path, output_dir);
    waitbar(i/n, f, sprintf('Batch Processing TAVA tEGG files... Progress: %d %%', floor(i/n*100)));
end
close(f)
