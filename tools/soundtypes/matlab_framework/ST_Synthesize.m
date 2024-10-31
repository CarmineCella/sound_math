function d_outWave = ST_Synthesize(d_dictionaryMat,ui_indexVec,st_pars,d_enVec)

[ui_atomLength,~,ui_numChannels] = size(d_dictionaryMat);
ui_numAtoms = length(ui_indexVec);

d_outWave = zeros(ui_atomLength*(ui_numAtoms+20)*st_pars.d_overlapFactor,ui_numChannels);  

if ~exist('d_enVec')
  d_enVec = ones(ui_numAtoms,1);
end
  

% overlap add
ui_first = 1;
for i=1:ui_numAtoms
    ui_last = ui_first + ui_atomLength - 1;
    d_outWave(ui_first:ui_last,:) = d_outWave(ui_first:ui_last,:) + squeeze(d_dictionaryMat(:,ui_indexVec(i),:)) * d_enVec(i);
    ui_first = ui_first + st_pars.ui_hopSize;
end



