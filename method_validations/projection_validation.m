
addpath('../projection/');


%% high value expected when it is completely in the spanned subspace
r = rand(3,11);
d = repmat(r, 30, 1);
d = d + 0.05 * randn(size(d));

bla = project(d, r, 1000, 0.5);
h   = plot_proj(bla);
figure(h);


%% small value expected for random data in random space
space   = rand(3,11);
dat     = rand(100, 11);
bla = project(dat, space, 1000, 0.5);
h   = plot_proj(bla);
figure(h);
