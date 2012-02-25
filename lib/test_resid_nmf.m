function [v var_res]=test_resid_nmf( matin , config)

% simple residual test

% this message is likely to happen in my case because of the highly
% structured data, espacially for the nonevoked data
warning off stats:nnmf:LowRank

centered = matin - repmat(mean(matin), size(matin, 1), 1);
SS = sum(sum(centered.^2));
v = NaN(config.Niter_exploration, min(size(matin)));
for i=1:min(size(matin))
    for j=1:config.Niter_exploration
        [w,h]=nnmf(matin',i);
        D = matin' - w*h;
        v(j,i) = sum(sum(D.^2));
    end
end
var_res = std(v);
v = (1 - (min(v)/SS)) * 100;

warning on stats:nnmf:LowRank
