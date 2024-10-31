function [d_transitionMat,str_NGrams,ui_ngramIndex] = ST_EstimateRules(ui_indexVec,b_doPlots,ui_markovOrder)

% markov chain
if ui_markovOrder==1
    d_transitionMat = hmmestimate(ui_indexVec,ui_indexVec);
    str_NGrams = [];
    ui_ngramIndex = [];
    N = max(ui_indexVec);
else
  
    N = max(ui_indexVec);
    
    % extract contiguous sequences of 2 items from the above
    for i=1:ui_markovOrder
        n_gramMat(i,:) = ui_indexVec(i:end-(ui_markovOrder-(i-1)));
    end
    str_NGrams = num2str( n_gramMat' );  % !
    ui_ngramIndex = grp2idx(str_NGrams);
    
    ui_numNGrams = length(unique(ui_ngramIndex));
    
    % items following the bigrams
    ui_follow = ui_indexVec(ui_markovOrder+1:end);
    
    d_transitionMat = full(sparse(ui_ngramIndex,ui_follow,1,ui_numNGrams,N));
    
    % normalize (unity sum)
    d_transitionMat = d_transitionMat./repmat(sum(d_transitionMat,2),1,N);
end


if b_doPlots
    figure(2);
    imagesc(d_transitionMat);
end