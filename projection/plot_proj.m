function h = plot_proj(in)

% a function to plot the results of a projection anaylisis and labels it.
% labelling has to be changed when for example test data is switched on in
% './project.m'

h = figure('Visible', 'off');

bins = 100;

ratios      = NaN(4,length(in.shuf_ratios));
ratios(1,:) = in.rand_ratios;
ratios(2,:) = in.shuf_ratios;
ratios(3,:) = in.test_ratios;
ratios(4,:) = in.rcov_ratios;

text = {'random data (same mean)', 'shuffled data', 'random data (same mean in channels)', 'shuffled data (only in channels)'};

for i = 1:size(ratios,1)
    subplot(4,1,i);
    hist(ratios(i,:),bins);
    n = hist(ratios(i,:),bins);
    title(text{i});
    xlim([0 max([mean(in.test_ratios) in.data_ratio]) * 1.2]);
    hold on
    plot(in.data_ratio * ones(1,max(n)), 1:max(n), '*r');
    hold off;
end


