function ui_indexGenVec = ST_GenerateSeq(d_transitionMat,ui_numFrames,ui_markovOrder,str_NGrams,ui_ngramIndex)


if ui_markovOrder==1
  

  % first non-zero-prob row
  ui_keepInd = find(sum(d_transitionMat,2)~=0);
  if ui_keepInd(1)>1
    d_transitionMat(1,ui_keepInd(1)) = 1;
  end
  
  ui_indexGenVec = hmmgenerate(ui_numFrames,d_transitionMat,eye(size(d_transitionMat)));
  
else
  
  % cumulative probabilities
  d_transitionMat = cumsum(d_transitionMat,2);
  numStates = size(d_transitionMat,2);
  
  statechange = rand(1,ui_numFrames);
  
  str_NGrams = str2num(str_NGrams);
  
  currentstate = 1;
  currNgram = str_NGrams(currentstate,:);
  
  for i=1:ui_numFrames
    
    stateVal = statechange(i);
    state = 1;
    for nextSymbol = numStates-1:-1:1
        if stateVal > d_transitionMat(currentstate,nextSymbol)
            state = nextSymbol + 1;
            break;
        end
    end
    
    ui_indexGenVec(i) = state;
    
    % map to n-grams
    str_nextNGram = [currNgram(2:end) state];  
    for j=1:size(str_NGrams,1)
      if all(str_NGrams(j,:) == str_nextNGram)
        break;
      end
    end
    
    if j==size(str_NGrams,1)   % n-gram not found
      currNgram = str_NGrams(1,:);  % restart from first
      currentstate = ui_ngramIndex(1);
    else
      currNgram = str_nextNGram;
      currentstate = ui_ngramIndex(j);
    end
    
  end
  
end