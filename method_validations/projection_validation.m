
% to test the validity of our projection method

r = rand(3,11);
d = repmat(r, 30, 1);
d = d + 0.05 * randn(size(d));

bla = project(d, r, 1000, 0.5);
h = plot_proj(bla);
figure(h);
