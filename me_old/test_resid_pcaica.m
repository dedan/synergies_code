function v=test_resid_pcaica( matin, config)

SS = sum(sum(matin));
m  = mean(matin);
lm = size(matin,1);

N = config.n_used_channels;

[w,h]=princomp(matin);
v = zeros(1,N);

for i=1:N
   dif = matin - ones(lm,1)*m - (w(:,1:i)*h(:,1:i)')';
   v(i) = sum(sum(abs(dif)));
end
v = (v/SS) * 100;
