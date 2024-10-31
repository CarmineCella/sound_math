function [d_stft,st_pars,d_atoms] = ST_STFT(d_waveForm,st_pars)

st_pars.ui_winSize = 2^nextpow2(st_pars.d_fs*st_pars.d_winLength);
st_pars.ui_hopSize = round(st_pars.ui_winSize*st_pars.d_overlapFactor);

ui_numFrames = floor(size(d_waveForm,1)/st_pars.ui_hopSize) - 1;

ui_numChannels = size(d_waveForm,2);

d_atoms = zeros(st_pars.ui_winSize,ui_numFrames,ui_numChannels);

d_winVec = repmat(blackman(st_pars.ui_winSize),1,ui_numChannels);

ui_pos = 1;
for fi=1:ui_numFrames - 3
    d_atoms(:,fi,:) = d_waveForm(ui_pos:ui_pos+st_pars.ui_winSize-1,:).*d_winVec;  
    ui_pos = ui_pos + st_pars.ui_hopSize - 1;
end

% force mono for fft
d_atomsFFT = mean(d_atoms,3);

d_stft = fft(d_atomsFFT,st_pars.ui_winSize);

% actually, spectrogram
d_stft(round(st_pars.ui_winSize/2)+1:end,:) = [];
d_stft = abs(d_stft);

