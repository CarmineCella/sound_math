function d_dictionaryMat = ST_CreateTypes(d_atoms,ui_indexVec,ui_abstraction,d_distVec)

ui_atomLength = size(d_atoms,1);
d_dictionaryMat = zeros(ui_atomLength,ui_abstraction,size(d_atoms,3));
d_distVec = exp(-d_distVec');

% weighted sum method
for i=1:ui_abstraction
    d_weights = repmat(d_distVec(ui_indexVec==i),ui_atomLength,1);
    for ci=1:size(d_atoms,3)
      d_dictionaryMat(:,i,ci) = sum(d_atoms(:,ui_indexVec==i,ci).*d_weights,2);
    end
end
