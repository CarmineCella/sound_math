% main script for sound types analysis and synthesis

close all;

% SOURCE SIGNAL
% ====================================================================

% str_sourceSound = './input_sounds/Beethoven_Symph7.wav';
% str_sourceSound = './input_sounds/beet_piano1.wav';
% str_sourceSound = './input_sounds/rien.wav';
% str_sourceSound = './input_sounds/guitar2.wav';
% str_sourceSound = './input_sounds/trumpet.wav';
% str_sourceSound = './input_sounds/TenorSax.wav';
% str_sourceSound = './input_sounds/cage.wav';
% str_sourceSound = './input_sounds/brad.wav';
% str_sourceSound = './input_sounds/APM1.wav';
% str_sourceSound = './input_sounds/piano_drums.wav';
% str_sourceSound = './input_sounds/god.wav';
str_sourceSound = '../../../samples/voices/god_vocal.wav';
% str_sourceSound = './input_sounds/nude_voice.wav';



str_outPath = './output_sounds';
[~,str_sourceOutName] = fileparts(str_sourceSound);
str_testOut = fullfile(str_outPath,[str_sourceOutName,'_out_SOURCE']);

% analysis parameters
st_pars.d_winLength = .16;
st_pars.d_overlapFactor = 1/4;
st_pars.d_clusterFactor = .9;   % number of clusters as percentage of frames
st_pars.ui_markovOrder = 1;

st_pars.b_doPlots = 0;

% PROCESS SOURCE SIGNAL
[d_sourceDictMat,d_sourceCentroidMat,~,st_pars,st_normPars.source,d_sourceTransitionMat,str_sourceNGrams,ui_sourceNgramIndex] = ST_process(str_sourceSound,st_pars,str_testOut);

% TARGET SIGNAL
% ====================================================================

% str_targetSound = './input_sounds/beet_piano1.wav';
% str_targetSound = './input_sounds/beet_piano1_TEST.wav';
str_targetSound = '../../../samples/voices/cage.wav';
% str_targetSound = './input_sounds/written.wav';

[~,str_targetOutName] = fileparts(str_targetSound);
str_testOut = fullfile(str_outPath,[str_targetOutName,'_out_TARGET']);

% analysis parameters
st_pars.d_winLength = .16;
st_pars.d_overlapFactor = 1/4;
st_pars.d_clusterFactor = 0.9;   % number of clusters as percentage of frames
st_pars.ui_markovOrder = 1;

% PROCESS TARGET SIGNAL
[~,d_targetCentroidMat,ui_targetIndexVec,~,st_normPars.target,d_targetTransitionMat,~,~,d_targetEnVec] = ST_process(str_targetSound,st_pars,str_testOut);

% PERFORM MATCHING
% ====================================================================

st_pars.b_doPlots = 0;

str_testOut = fullfile(str_outPath,[str_sourceOutName,'_',str_targetOutName,'_out']);
d_distMat = ST_match(d_sourceDictMat,d_sourceCentroidMat,d_targetCentroidMat,ui_targetIndexVec,st_pars,str_testOut,st_normPars,d_targetEnVec);

% PROBABILITY MERGING
% ====================================================================
ui_numTargetFrames = length(ui_targetIndexVec);
ST_probMerge(d_sourceDictMat,d_sourceTransitionMat,d_targetTransitionMat,d_distMat,st_pars,str_sourceNGrams,ui_sourceNgramIndex,ui_numTargetFrames,str_testOut);
