function ST_probMerge(d_sourceDictMat,d_sourceTransitionMat,d_targetTransitionMat,d_distMat,st_pars,str_sourceNGrams,ui_sourceNgramIndex,ui_numTargetFrames,str_testOut)

d_mergeFactor = .9; % 1: 100% target

% tables of closest centroids

% source to target
% [~,ui_S2TcentroidMap] = min(d_distMat,[],2);

% target to source
% [~,ui_T2ScentroidMap] = min(d_distMat,[],1);

if size(d_targetTransitionMat,1) < size(d_sourceTransitionMat,1)
  ui_T2ScentroidMap = [];
  for i=1:size(d_targetTransitionMat,1)
    [~,curr_indVec] = sort(d_distMat(:,i));
    for j=1:length(curr_indVec)
      if isempty(find(curr_indVec(j) == ui_T2ScentroidMap))
        ui_T2ScentroidMap = [ui_T2ScentroidMap curr_indVec(j)];
        break;
      end
    end
  end
  
end

% target to source merging
d_sourceTransitionMERGEDMat = zeros(size(d_sourceTransitionMat));
d_sourceTransitionMERGEDMat(ui_T2ScentroidMap,ui_T2ScentroidMap) = d_targetTransitionMat;

d_sourceTransitionMERGEDMat = d_sourceTransitionMERGEDMat*d_mergeFactor + d_sourceTransitionMat*(1-d_mergeFactor);

% normalize (unity sum)
d_sourceTransitionMERGEDMat = d_sourceTransitionMERGEDMat./repmat(sum(d_sourceTransitionMERGEDMat,2),1,size(d_sourceTransitionMERGEDMat,1));

d_sourceTransitionMERGEDMat(isnan(d_sourceTransitionMERGEDMat)) = 0;

ui_indexGenVec = ST_GenerateSeq(d_sourceTransitionMERGEDMat,ui_numTargetFrames,st_pars.ui_markovOrder,str_sourceNGrams,ui_sourceNgramIndex);
d_outWave = ST_Synthesize(d_sourceDictMat,ui_indexGenVec,st_pars);

d_outWave = d_outWave/max(max(abs(d_outWave)))*.999;

wavwrite(d_outWave,st_pars.d_fs,st_pars.ui_bits,[str_testOut,'_PROBMERGED.wav']);