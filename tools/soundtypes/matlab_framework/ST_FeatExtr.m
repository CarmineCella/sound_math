function [d_featMat,st_normPars,d_enVec] = ST_FeatExtr(d_spec,st_pars)

ui_numFrames = size(d_spec,2);

ui_numMFCC = 12;

st_feat.centroid.d_value = zeros(ui_numFrames,1);
st_feat.skewness.d_value   = zeros(ui_numFrames,1);
st_feat.mfcc.d_value = zeros(ui_numFrames,ui_numMFCC);

d_freqVec = st_pars.d_fs/st_pars.ui_winSize*(1:round(st_pars.ui_winSize/2))';

% extract features
for fi=1:ui_numFrames
  % centroid
  st_feat.centroid.d_value(fi) = sum(d_freqVec.*d_spec(:,fi))/sum(d_spec(:,fi));
  
  % flatness
  st_feat.flatness.d_value(fi) = geomean(d_spec(:,fi))/mean(d_spec(:,fi));
  
  % skewness
  st_feat.skewness.d_value(fi) = skewness(d_spec(:,fi));
  
end

% mfcc
st_feat.mfcc.d_value = ST_mfcc(d_spec,st_pars.d_fs,ui_numMFCC);

% create feature matrix
d_featMat(1,:) = [st_feat.centroid.d_value];
d_featMat(2,:) = [st_feat.flatness.d_value];
d_featMat(3,:) = [st_feat.skewness.d_value];
d_featMat(4:4+ui_numMFCC-1,:) = st_feat.mfcc.d_value;

ui_numFeat = size(d_featMat,1);

% % weight features by energy (?)
d_enVec = sum(d_spec,1);
% d_featMat = d_featMat.*repmat(d_enVec,ui_numFeat,1);

% normalize features
st_normPars.mean = nanmean(d_featMat,2);
st_normPars.std = nanstd(d_featMat,1,2);
d_featMat = (d_featMat - repmat(st_normPars.mean,1,ui_numFrames))./repmat(st_normPars.std,1,ui_numFrames);

