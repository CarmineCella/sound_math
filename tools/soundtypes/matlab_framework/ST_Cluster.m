function [ui_indexVec,d_centroidMat,d_distVec,ui_elimIndex] = ST_Cluster(d_featMat,ui_numTypes,b_doPlots)

% using kmeans
[ui_indexVec,d_centroidMat,~,d_distVec] = kmeans(d_featMat',ui_numTypes,'emptyaction','drop');

% vector of distance of each point to its centroid
d_distVec = min(d_distVec,[],2);

ui_elimIndex = find(isnan(ui_indexVec));
ui_indexVec(ui_elimIndex) = [];

if b_doPlots
    str_colors = 'rgbcmyk';
    figure(1);
    for i=1:ui_numTypes
        scatter3(d_featMat(1,ui_indexVec==i),d_featMat(2,ui_indexVec==i),d_featMat(3,ui_indexVec==i),[str_colors(mod(i-1,7)+1),'.']);
        hold on;
        scatter3(d_centroidMat(:,1),d_centroidMat(:,2),d_centroidMat(:,3),50,...
            'filled','LineWidth',2,'MarkerFaceColor',str_colors(mod(i-1,7)+1),'MarkerEdgeColor','k');
    end
end