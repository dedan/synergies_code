function h = plot_proj(in)

% a function to plot the results of a projection anaylisis and labels it.
% labelling has to be changed when for example test data is switched on in
% './project.m'

h = figure('Visible', 'off');

bins = 100;

ratios      = NaN(4,length(in.shuf_ratios));
ratios(1,:) = in.rand_ratios;
ratios(2,:) = in.shuf_ratios;
ratios(3,:) = in.rand_chan_ratios;
ratios(4,:) = in.shuf_chan_ratios;

text = {'random data (same mean)', 'shuffled data', 'random data (same mean in channels)', 'shuffled data (only in channels)'};

for i = 1:size(ratios,1)
    subplot(5,1,i+1);
    title(text{i});
    hold on
    [n, xout] = hist(ratios(i,:),bins);
    bar(xout, n / size(ratios,2)  , 'b');
    plot(mean(in.ratio_dist) * ones(1,10), linspace(0,max(n / size(ratios,2)), 10), '*r');
    hold off;
    xlim([0 max([mean(ratios(1,:)) mean(in.ratio_dist)]) + 2 * std(in.ratio_dist)]);
end

subplot(5,1,1);
[n,xout] = hist(in.ratio_dist, bins);
bar(xout, n / length(in.ratio_dist), 'r');
xlim([0 max([mean(ratios(1,:)) mean(in.ratio_dist)]) + 2 * std(in.ratio_dist)]);
title('distribution of projection ratios')