function d_mfccMat = ST_mfcc(d_spec,d_Fs,ui_numCoeff)

ui_numFilt = 40;

% compute PSD
d_spec = d_spec.^2;

% create filter bank
ui_halfFFTlen = size(d_spec,1);                                          
d_HzPerBin = d_Fs/(2*ui_halfFFTlen);                                             
d_minMel = freq2mel(d_HzPerBin);                                   
d_maxMel = freq2mel(ui_halfFFTlen*d_HzPerBin);                                   
d_melPerBin = (d_maxMel-d_minMel)/(ui_numFilt+1);                                 

d_centerMelVec = [d_minMel,d_melPerBin*ones(1,ui_numFilt+1)];
d_centerMelVec = round(mel2freq(cumsum(d_centerMelVec))/d_HzPerBin);

% create triangular filters
d_filtBankMat = zeros(ui_numFilt,ui_halfFFTlen);
for i = 1:ui_numFilt
  d_filtBankMat(i,d_centerMelVec(i):d_centerMelVec(i+1))   = linspace(0,1,d_centerMelVec(i+1)-d_centerMelVec(i)+1);
  d_filtBankMat(i,d_centerMelVec(i+1):d_centerMelVec(i+2)) = linspace(1,0,d_centerMelVec(i+2)-d_centerMelVec(i+1)+1);
end

d_mfccMat = dct(log(d_filtBankMat*d_spec+eps));

% ignore energy
d_mfccMat = d_mfccMat(2:ui_numCoeff+1,:);


function d_mel = freq2mel(d_freq)

d_mel = 1127.01048*log(1+d_freq/700);

function d_freq = mel2freq(d_mel)

d_freq = 700*(exp(d_mel/1127.01048)-1);