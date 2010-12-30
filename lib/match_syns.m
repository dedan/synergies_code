%This function matches rows from W1 and W2 and sorts them according 

%scores gives the correlation coefficient between columns of W1 and their
%matched columns in W2. 

function [w1, w2, scores] = match_syns(w1,w2, baseline)



s       = size(w1,1);
scores  = zeros(1,s);
idx     = zeros(1,s);

% similarity of synergies is measured by dot product
w1      = normr(w1);
w2      = normr(w2);
dot_abs = abs(w1*w2');
dots    = (w1*w2');


for i=1:s
    [m,irow]     = max(dot_abs);
    [~,icol]     = max(m);
    irow         = irow(icol);

    % turn first syn if dot product negative
    if dots(irow, icol) < 0
        w1(irow,:) = w1(irow,:) * -1;
    end

    dot_abs(irow,:) = -Inf;
    dot_abs(:,icol) = -Inf;
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
if nargin == 3
    scores = (scores - baseline) / (1 - baseline);
end 



