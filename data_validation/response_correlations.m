
% as many times before, the issue or normalization came up.
% We discussed whether we should normalize the data from the
% natural movements, if yes -> how and also whether to normalize the evoked
% responses. Can we do this? Are the differences between the channels (some of them
% are more active then others) correlated for natural movement data and evoked responses.
% this file computes the correlation between the two and as a result showed
% that means and stds are not correlated.


data_path = '~/projects/yifat_paper/results/data/';
monks = {'chalva', 'vega'};

for i = 1:length(monks)
    load(['~/projects/yifat_paper/results/data_without_norm/' 'all_data_' monks{i}]);
    load([data_path 'nat_mov_res_' monks{i}])
    c2take = nat_mov_res.c2take;

    all_natural = [];
    for j = 1:length(sessions)
        all_natural = [all_natural; sessions(j).mats(1).data_raw];
    end

    means_nat = mean(all_natural);
    stds_nat = std(all_natural);

    load([data_path 'evoked_data_' monks{i}]);
    all_resps = vertcat(resps.response);
    means_evo = mean(all_resps);
    stds_evo = std(all_resps);

    figure(i)
    subplot 211
    plot(means_nat(c2take), means_evo, '.');
    [r, p] = corrcoef(means_nat(c2take), means_evo)
    title(['means ' monks{i}])

    subplot 212
    plot(stds_nat(c2take), stds_evo, '.');
    [r, p] = corrcoef(stds_nat(c2take), stds_evo)
    title(['stds ' monks{i}])

end