% test of variability maximitaion in sound types matching
close all;

ui_numPoints = 40;
SdB = mvnrnd([2 2],[.05 0 ; 0 .05],ui_numPoints); % source database
TdB = mvnrnd([-2 -2],[1 .5 ; .5 1],ui_numPoints); % target database

% compute means and stds
mu_SdB = mean(SdB);
sigma_SdB = std(SdB);
mu_TdB = mean(TdB);
sigma_TdB = std(TdB);
mu_total = mean([SdB ; TdB]);
sigma_total = std([SdB ; TdB]);

figure(1);
subplot(321); axis tight;
scatter(SdB(:,1),SdB(:,2),'ok');
hold on; box on; grid on; 
scatter(TdB(:,1),TdB(:,2),'+r');
title('1. original data');

% separate normalization
SdB_norm1 = (SdB - repmat(mu_SdB,ui_numPoints,1))./repmat(sigma_SdB,ui_numPoints,1);
TdB_norm1 = (TdB - repmat(mu_TdB,ui_numPoints,1))./repmat(sigma_TdB,ui_numPoints,1);

subplot(322);
scatter(SdB_norm1(:,1),SdB_norm1(:,2),'ok');
hold on; box on; grid on;
scatter(TdB_norm1(:,1),TdB_norm1(:,2),'+r');
title('2. separate normalization');

% cross normalization
SdB_norm2 = (SdB - repmat(mu_TdB,ui_numPoints,1))./repmat(sigma_TdB,ui_numPoints,1);
TdB_norm2 = (TdB - repmat(mu_TdB,ui_numPoints,1))./repmat(sigma_TdB,ui_numPoints,1);

subplot(323);
scatter(SdB_norm2(:,1),SdB_norm2(:,2),'ok');
hold on; box on; grid on;
scatter(TdB_norm2(:,1),TdB_norm2(:,2),'+r');
title('3. cross normalization');

% separate normalization putting back the means
SdB_norm3 = (SdB - repmat(mu_SdB,ui_numPoints,1))./repmat(sigma_SdB,ui_numPoints,1) + repmat(mu_SdB,ui_numPoints,1);
TdB_norm3 = (TdB - repmat(mu_TdB,ui_numPoints,1))./repmat(sigma_TdB,ui_numPoints,1) + repmat(mu_TdB,ui_numPoints,1);

subplot(324);
scatter(SdB_norm3(:,1),SdB_norm3(:,2),'ok');
hold on; box on; grid on;
scatter(TdB_norm3(:,1),TdB_norm3(:,2),'+r');
title('4. separate normalization putting back means');

% separate normalization with s->t cross-renormalization
SdB_norm4 = (SdB - repmat(mu_SdB,ui_numPoints,1))./repmat(sigma_SdB,ui_numPoints,1);
TdB_norm4 = (TdB - repmat(mu_TdB,ui_numPoints,1))./repmat(sigma_TdB,ui_numPoints,1);
TdB_norm4 = (TdB_norm4 - repmat(mu_SdB,ui_numPoints,1))./repmat(sigma_SdB,ui_numPoints,1);

subplot(325);
scatter(SdB_norm4(:,1),SdB_norm4(:,2),'ok');
hold on; box on; grid on;
scatter(TdB_norm4(:,1),TdB_norm4(:,2),'+r');
title('5. separate normalization with s->t cross-renormalization');

% separate normalization with t->s cross-renormalization
SdB_norm5 = (SdB - repmat(mu_SdB,ui_numPoints,1))./repmat(sigma_SdB,ui_numPoints,1);
TdB_norm5 = (TdB - repmat(mu_TdB,ui_numPoints,1))./repmat(sigma_TdB,ui_numPoints,1);
SdB_norm5 = (SdB_norm5 - repmat(mu_TdB,ui_numPoints,1))./repmat(sigma_TdB,ui_numPoints,1);

subplot(326);
scatter(SdB_norm5(:,1),SdB_norm5(:,2),'ok');
hold on; box on; grid on;
scatter(TdB_norm5(:,1),TdB_norm5(:,2),'+r');
title('6. separate normalization with t->s cross-renormalization');