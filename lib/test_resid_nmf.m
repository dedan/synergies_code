function [v var_res]=test_resid_nmf( matin , config)

% simple residual test

SS = sum(sum(matin));
v = NaN(config.Niter_exploration, min(size(matin)));
for i=1:min(size(matin))
    for j=1:config.Niter_exploration
        [w,h]=nnmf(matin',i);
        D = matin' - w*h;
        v(j,i) = sum(sum(abs(D)));
    end
end
var_res = std(v);
v = (min(v)/SS) * 100;
    