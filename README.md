Code for the method described in our 2023 ICASSP paper: https://ccrma.stanford.edu/~cnoufi/demos/TAVA-PhonemeRemoval/

# TAVA-LanguagePerceptionRemoval

### How to use:

The wrapper script is:
main/generate_stimulus.m 

To generate a stimulus from a single audio file, run from Matlab's command line:
generate_stimulus(‘complete/path/to/audiofile.wav’)

The input audio file must have the channel structure (mono_speech.wav, EGG.wav, ...)
