% This function matches rows from W1 and W2 and sorts them according

% scores gives the correlation coefficient between rows of W1 and their
% matched rows in W2.  
% type:  1 dot product matching  -- 2 corrcoef matching

function [w1, w2, scores, idx] = match_syns(w1,w2, type, baseline)

s       = size(w1,1);
scores  = zeros(1,s);
idx     = zeros(1,s);
w1      = normr(w1);
w2      = normr(w2);


% similarity of synergies is measured by dot product
if type == 1

    dist_abs    = abs(w1*w2');
    dist        = (w1*w2');
    
% measured by corrcoef
else
    
    c           = corrcoef([w1' w2']);
    % omit the comparisons with itself
    c           = triu(c(1:s,s+1:end));
    c           = c + tril(c',-1);
    dist_abs    = abs(c);
    dist        = c;
end

% do the matching      
for i=1:s
    [m,irow]     = max(dist_abs);
    [~,icol]     = max(m);
    irow         = irow(icol);
    
    % turn first syn if dot product negative
    if dist(irow, icol) < 0
        w1(irow,:) = w1(irow,:) * -1;
    end
    
    dist_abs(irow,:) = -Inf;
    dist_abs(:,icol) = -Inf;
    idx(irow)       = icol;
end

% sort the vectors in w2 according to w1
w2 = w2(idx,:);

% compute scores by correlation coefficients
for i = 1:s
    c           = corrcoef(w1(i,:), w2(i,:));
    scores(i)   = c(2,1);
end

% normalize scores with estimated baseline
if nargin == 4
    scores = (scores - baseline) / (1 - baseline);
end



