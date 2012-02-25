function v=test_resid_pcaica( matin, config)

% simple residual test

centered = matin - repmat(mean(matin), size(matin, 1), 1);
SS = sum(sum(centered.^2));
m  = mean(matin);
lm = size(matin,1);

[w,h]=princomp(matin);
v = zeros(1,min(size(matin)));

for i=1:min(size(matin))
    % mean has to be used for normalization, is removed during princomp
   dif = matin - ones(lm,1)*m - (w(:,1:i)*h(:,1:i)')';
   v(i) = sum(sum(dif.^2));
end
v = (v/SS) * 100;
