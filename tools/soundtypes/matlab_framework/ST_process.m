function [d_dictionaryMat,d_centroidMat,ui_indexVec,st_pars,st_normPars,d_transitionMat,str_NGrams,ui_ngramIndex,d_enVec] = ST_process(str_soundFile,st_pars,str_testOut)


% read sound file
[d_waveForm,st_pars.d_fs,st_pars.ui_bits] = wavread(str_soundFile);

% normalize
d_waveForm = d_waveForm/max(max(abs(d_waveForm)));

% compute spectrogram
[d_spec,st_pars,d_atoms] = ST_STFT(d_waveForm,st_pars);

% feature extraction
[d_featMat,st_normPars,d_enVec] = ST_FeatExtr(d_spec,st_pars);

% cluster
ui_numFrames = size(d_featMat,2);
st_pars.ui_numClusters = round(st_pars.d_clusterFactor*ui_numFrames);
[ui_indexVec,d_centroidMat,d_distVec,ui_elimIndex] = ST_Cluster(d_featMat,st_pars.ui_numClusters,st_pars.b_doPlots);
d_atoms(:,ui_elimIndex,:) = [];

% estimate transition probabilities
[d_transitionMat,str_NGrams,ui_ngramIndex] = ST_EstimateRules(ui_indexVec,st_pars.b_doPlots,st_pars.ui_markovOrder);

% create sound types from clusters
d_dictionaryMat = ST_CreateTypes(d_atoms,ui_indexVec,st_pars.ui_numClusters,d_distVec);

% represent input signal with clusters (REBUILD)
d_outWave = ST_Synthesize(d_dictionaryMat,ui_indexVec,st_pars);

d_outWave = d_outWave/max(max(abs(d_outWave)))*.999;
wavwrite(d_outWave,st_pars.d_fs,st_pars.ui_bits,[str_testOut,'_REBUILD.wav']);

% GENERATE
ui_indexGenVec = ST_GenerateSeq(d_transitionMat,ui_numFrames,st_pars.ui_markovOrder,str_NGrams,ui_ngramIndex);
d_outWave = ST_Synthesize(d_dictionaryMat,ui_indexGenVec,st_pars);

d_outWave = d_outWave/max(max(abs(d_outWave)))*.999;

wavwrite(d_outWave,st_pars.d_fs,st_pars.ui_bits,[str_testOut,'_GENERATE.wav']);