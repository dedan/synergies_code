
set_size = 1000;

orig = rand(4,12);

data = zeros(set_size,12);

for i = 1:set_size
   data(i,:) = rand(1,4) * orig;    
end

data = data + 0.4 *randn(size(data));

imagesc(data);

%[W H] = nnmf(data,4);

for i=1:7
    [lambda,psi,T,stats] = factoran(data,i, 'scores', 'wls');
    res(i) = stats.loglike;
end
plot(res);

