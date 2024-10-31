function d_distMat = ST_match(d_sourceDictMat,d_sourceCentroidMat,d_targetCentroidMat,ui_targetIndexVec,st_pars,str_testOut,...
  st_normPars,d_targetEnVec)

% str_normMethod = 'none';
str_normMethod = 'separate';
% str_normMethod = 'separate_mean_back';  % equivalent to not removing mean (normalize only with std)

ui_numSourceCentroids = size(d_sourceCentroidMat,1);
ui_numTargetCentroids = size(d_targetCentroidMat,1);

% apply normalizations
switch str_normMethod
  case 'none'
    % de-normalize
    d_sourceCentroidMat = d_sourceCentroidMat.*repmat(st_normPars.source.std',ui_numSourceCentroids,1) + ...
      repmat(st_normPars.source.mean',ui_numSourceCentroids,1);
    d_targetCentroidMat = d_targetCentroidMat.*repmat(st_normPars.target.std',ui_numTargetCentroids,1) + ...
      repmat(st_normPars.target.mean',ui_numTargetCentroids,1);
  case 'separate'
    % do nothing: already normalized separately (for clustering)
  case 'separate_mean_back'
    d_sourceCentroidMat = d_sourceCentroidMat + repmat(st_normPars.source.mean',ui_numSourceCentroids,1)./...
      repmat(st_normPars.source.std',ui_numSourceCentroids,1);
    d_targetCentroidMat = d_targetCentroidMat + repmat(st_normPars.target.mean',ui_numTargetCentroids,1)./...
      repmat(st_normPars.target.std',ui_numTargetCentroids,1);
end



if st_pars.b_doPlots
  scatter3(d_sourceCentroidMat(:,1),d_sourceCentroidMat(:,2),d_sourceCentroidMat(:,3),'k.');
  hold on;
  scatter3(d_targetCentroidMat(:,1),d_targetCentroidMat(:,2),d_targetCentroidMat(:,3),'r.');
  legend({'source','target'});
end

% create distance matrix

d_distMat = zeros(ui_numSourceCentroids,ui_numTargetCentroids);

for si=1:ui_numSourceCentroids
  for ti=1:ui_numTargetCentroids
    d_distMat(si,ti) = norm(d_sourceCentroidMat(si,:) -  d_targetCentroidMat(ti,:)); % euclidean
  end
end

% select centroids from the source closest to the target sequence
for i=1:length(ui_targetIndexVec)
  [~,ui_sourceIndexVec(i)] = min(d_distMat(:,ui_targetIndexVec(i)));
end

% plot trajectories
if 0
  for i=1:length(ui_targetIndexVec)
    scatter3(d_sourceCentroidMat(ui_sourceIndexVec(i),1),d_sourceCentroidMat(ui_sourceIndexVec(i),2),d_sourceCentroidMat(ui_sourceIndexVec(i),3),'go','LineWidth',6);
    scatter3(d_targetCentroidMat(ui_targetIndexVec(i),1),d_targetCentroidMat(ui_targetIndexVec(i),2),d_targetCentroidMat(ui_targetIndexVec(i),3),'mo','LineWidth',6);
    pause
  end
end

% replace target dictionary atoms with source atoms
d_outWave = ST_Synthesize(d_sourceDictMat,ui_sourceIndexVec,st_pars,d_targetEnVec);

d_outWave = d_outWave/max(max(abs(d_outWave)));

wavwrite(d_outWave,st_pars.d_fs,st_pars.ui_bits,[str_testOut,'_MATCHED.wav']);