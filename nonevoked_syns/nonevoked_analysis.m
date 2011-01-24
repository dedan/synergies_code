
%% file for the analysis of natural movement

% it is basically the same as the nonevoked syns script, but now for all
% three monkey, some of them don't have both handpositions, and so on. It
% should in the end produce all the data (concerning natural movement)
% which might be needed for an article or poster.

% procedures which show the validity of methods (like nmf stability) that
% do not necessarily rely on all data and the latest data will move to
% nonevoked_sessions

% for the beginning I will structure this file with analysis like done in
% the poster for Ein Gedi

%% configurations

clear;
clc;
conf = struct;
conf.opt                = statset('MaxIter',5);
conf.outpath            = '~/Documents/uni/yifat_lab/results/natural_mov/';
conf.inpath             = '/Volumes/LAB/results/data/';   
conf.names              = {'chalva', 'vega', 'darma'};
conf.n_monks            = length(conf.names);
conf.significant        = 25;
conf.max_channels       = 16;
conf.dim                = 3;
conf.Niter_exploration  = 2;
conf.n_best             = 2;
conf.n_trials           = 50;
conf.norm_matchscore    = false;
conf.image_format       = 'jpg';
conf.n_pd_unstable      = 4;
conf.only_sig_pd        = true;
conf.n_bootstrap        = 100;
mymap = [linspace(145/255,178/255,32)' linspace(167/255,213/255,32)' linspace(216/255,111/255,32)'];
mymap = [mymap; linspace(178/255,251/255,32)' linspace(213/255,147/255,32)' linspace(111/255,24/255,32)'];
conf.map = mymap;

delete([conf.outpath 'log.txt']);
diary([conf.outpath 'log.txt']);
clear mymap

res = struct;


%% load data
sessions = struct([]);
for i = 1:conf.n_monks
    if isempty(sessions)
        load([conf.inpath 'all_data_' conf.names{i}]);
%        res.(conf.names{i}).stats = stats;
    else
        tmp         = load([conf.inpath 'all_data_' conf.names{i}]);
        sessions    = [sessions tmp.sessions]; %#ok<AGROW>
%        res.(conf.names{i}).stats = tmp.stats;        
    end        
end

addpath('../lib'); 
addpath('../plots');
clear tmp stats i 



%% statistics
% First lets get some statistics on the Data we have available

for i = 1:length(sessions)
    sessions(i).target1     = size(sessions(i).mats(1).data,1); 
    sessions(i).n_channels  = length(find(sessions(i).channels)); 
    
    if length(sessions(i).mats) == 2
        sessions(i).target2    = size(sessions(i).mats(2).data,1); 
    else
        sessions(i).target2    = -1; 
    end
end

h = figure('Visible', 'off');
disp(['data from in total: ' num2str(length(sessions)) ' sessions']);

for i = 1:conf.n_monks
    idx.(conf.names{i})    = strcmp(conf.names{i}, {sessions.monk});
    
    disp([conf.names{i} ': ' num2str(length(find(idx.(conf.names{i})))) ' sessions' ]);

    % distribution of targets
    subplot(4,conf.n_monks,i);
    [y x] = hist([sessions(idx.(conf.names{i})).target1],-1:8);
    bar(x,y,'r');
    hold on
    [y x] = hist([sessions(idx.(conf.names{i})).target2],-1:8);
    bar(x+0.2,y,'b');
    hold off
    title(conf.names{i});
    xlabel('targets');

    % distribution of used channels
    subplot(4, conf.n_monks, conf.n_monks +i);
    [y x] = hist([sessions(idx.(conf.names{i})).n_channels],1:16);
    bar(x,y,'b');
    xlabel('channels');
    
    % distribution of trials
    subplot(4, conf.n_monks, 2*conf.n_monks +i);
    hist([sessions(idx.(conf.names{i})).trials], 100);
    xlabel('trials');
end
disp(' ');
disp(['mean number of trials: ' num2str(mean([sessions.trials]))]);

subplot(4,conf.n_monks, conf.n_monks*3+1);
hist([sessions.hands]);
title('handpositions');

saveas(h, [conf.outpath  'statistics.' conf.image_format]);
close(h);
clear y x h i 


%% sorting

for i = 1:conf.n_monks
    monk = conf.names{i};
    
    % sort out channels with less then 50 trials
    disp([monk ': sort out ' num2str(length(find([sessions(idx.(monk)).trials] <= conf.n_trials ))) ... 
        ' sessions because less then 50 trials']);
    idx.(monk) = idx.(monk) & [sessions.trials] > conf.n_trials;

    % sort out sessions with unstable pds
    % pronation pds because for all available
    data = zeros(length(find(idx.(monk))), conf.max_channels);
    inds = find(idx.(monk));
    for j = 1:length(inds)
        data(j,:) = sessions(inds(j)).pd(1,:); 
    end
    data(data > pi)  = data(data > pi) - 2* pi;
    
    all_mean = repmat(circ_mean(data), size(data,1), 1);
    all_std  = repmat(circ_std(data) * 2, size(data,1), 1);
    out      = sum(abs(data - all_mean) > all_std, 2)';
    
    disp([monk ': sort out ' num2str(length(find(out >= conf.n_pd_unstable))) ...
        ' sessions because unstable pds']);
    idx.(monk)(idx.(monk)) = idx.(monk)(idx.(monk)) & out < conf.n_pd_unstable;
end
clear data all_mean all_std i inds j monk out 



%% only use channels which are available for all sessions of a monk
for i = 1:conf.n_monks
    
    all_chan                     = vertcat(sessions(idx.(conf.names{i})).channels);
    res.(conf.names{i}).c2take   = all(all_chan);
end

clear all_chan i


%% stability of prefered directions (over sessions) (feather)
for i = 1:conf.n_monks

    max_hands   = max([sessions(idx.(conf.names{i})).hands]);
    ses2take    = idx.(conf.names{i}) & [sessions.hands] == max_hands;
    n2take      = length(find(ses2take));
    c2take_idx  = find(res.(conf.names{i}).c2take);
    n_feather   = ceil(length(c2take_idx)/2);
    
    sig     = 0.05;
    
    all_pd  = vertcat(sessions(ses2take).pd);
    all_p1  = vertcat(sessions(ses2take).p1);
    sigs    = all_p1 < sig;
    
    h = figure('Visible', 'off');
    
    if max_hands == 1
        ids(1).id = 1:n2take;
    else
        ids(1).id = 1:2:2*n2take-1;
        ids(2).id = 2:2:2*n2take;
    end
    
    % pd significance
    subplot(n_feather +2, 2, 1);
    agg = NaN(length(c2take_idx), max_hands);
    for j = 1:max_hands
        agg(:,j) = sum(all_p1(ids(j).id, c2take_idx) < 0.05)' ./ n2take * 100;
    end
    bar(c2take_idx,  agg);
    title('pd significance in percent of sessions');
    
    for j = 1:max_hands
        
        for k = 1:length(c2take_idx)
         
            tmp = all_pd(ids(j).id, :);
            [~, res.(conf.names{i}).cstd_all(j,k)] = circ_std(tmp(:, k));
            [~, res.(conf.names{i}).cstd_sig(j,k)] = circ_std(tmp(sigs(ids(j).id, k), k));                        
            res.(conf.names{i}).pds_all(j,k)       = circ_mean(tmp(:, k));            
            res.(conf.names{i}).pds_sig(j,k)       = circ_mean(tmp(sigs(ids(j).id, k), k));
        end
    end
    subplot(n_feather +2, 2, 2);
    bar(c2take_idx, res.(conf.names{i}).cstd_all(c2take_idx));
    title('cstds over all sessions');
    ylim([0 2]);
    
    subplot(n_feather +2, 2, 4);
    bar(c2take_idx, res.(conf.names{i}).cstd_sig(c2take_idx));
    title('cstds over significant sessions');
    ylim([0 2]);
    saveas(h, [conf.outpath  'pd_consist_bars_' conf.names{i} '.' conf.image_format]);
    close(h);    
    
    
    h = figure('Visible', 'off');    
    colors          = {'b', 'r'};
    for k = 1:length(c2take_idx)
        
        subplot(n_feather, 2, k);
        
        for j = 1:max_hands
            tmp   = all_pd(ids(j).id, c2take_idx(k));
            [x y] = pol2cart(tmp, ones(size(tmp)));
            hs = feather(x, y, colors{j});
            % overdraw the not significant pds with a dotted line
            set(hs(~sigs(ids(j).id, c2take_idx(k))), 'LineStyle', ':');            
            hold on
        end
        axis off
        title(['channel ' num2str(c2take_idx(k))]);
    end
        
    saveas(h, [conf.outpath  'pd_consist_feather_' conf.names{i} '.' conf.image_format]);
    close(h);
end



%% stability of prefered directions (over sessions) (scatter)
for i = 1:conf.n_monks

    h = figure('Visible', 'off');
    for j = 1:max_hands
        subplot(max_hands,1,j)
        inds = find(ses2take);
        for k = 1:length(inds)
            col = [k k k]/length(inds);
            in_deg = rad2deg(sessions(inds(k)).pd(j,c2take_idx));
            in_deg(in_deg < 0) = 360 + in_deg(in_deg < 0);
            plot(c2take_idx, in_deg, '^','MarkerSize',10, 'MarkerFaceColor', col);
            set(gca,'XTickLabel',arrayfun(@(x) sprintf('%.2f',x), ...
                res.(conf.names{i}).cstd_all(j,c2take_idx), 'UniformOutput', false))
            set(gca,'XTick', c2take_idx);
            hold on
        end
        grid
        xlabel('circular stds');
        ylabel('pds (deg)');
        hold off
    end
    saveas(h, [conf.outpath  'pd_consist_' conf.names{i} '.' conf.image_format]);
    close(h);
        
end


clear x y all_pd all_p1 all_p2 c2take_idx n_hands monk_first col in_deg inds id_pro 
clear id_sup bla h i j k max_hands n2take n_feather ses2take sig sigs tmp tmp1




%% rank1 analysis
% plot the rank1 values for all sessions and the different tests to see
% relation between the different tests. also check for difference according
% to handpos


% first store them in a convenient matrix
rank1 = NaN(4,length(sessions));
for i = 1:length(sessions)
    rank1(1,i)  = sessions(i).r_nmf_raw_pro(1);
    rank1(2,i)  = sessions(i).r_pca_raw_pro(1);
    rank1(3,i)  = sessions(i).r_nmf_raw_sup(1);
    rank1(4,i)  = sessions(i).r_pca_raw_sup(1);
end

h = figure('Visible', 'off');

% plot the rank1 values for pca and nmf and the line at the value at which 
% sessions are chosen as significant
subplot(2,2,1:2);
plot(rank1([1 3],:)','.');
hold on;
plot(ones(1,size(rank1, 2)) * conf.significant);

% separating lines between the monkeys
for j = 1:conf.n_monks
    plot(ones(1,10)*find(idx.(conf.names{j}), 1, 'last' ), 1:10:100, 'r--');
    text(find(idx.(conf.names{j}), 1 )+2, 70, conf.names{j});
end
hold off;
legend('nmf pro', 'nmf sup', 'Location', 'NorthWest');

% also calculate the correlation of the two rank1 values
subplot(2,2,3)
all_nmf = rank1([1 3],:);
all_pca = rank1([2 4],:);
[r p] = corrcoef(all_nmf(:), all_pca(:));
plot(all_nmf(:), all_pca(:), '.');
title(['nmf vs. pcaica, r: ' num2str(r(1,2)) ' - p: ' num2str(p(1,2))]);

subplot(2,2,4)
all_pro = rank1(1:2, rank1(4,:) > 0);
all_sup = rank1(3:4, rank1(4,:) > 0);
[r p] = corrcoef(all_pro(:), all_sup(:));
plot(all_pro(:), all_sup(:), '.');
title(['pro vs. sup, r: ' num2str(r(1,2)) ' - p: ' num2str(p(1,2))]);

saveas(h, [conf.outpath  'rank1.' conf.image_format]);    
close(h);


% sort out boring sessions
% we sort out sessions where the rank1 model already has an remaining error
% lower then 25 %
disp(' ')
disp('rank1 filtering');
for i = 1:conf.n_monks
    n_sort_out = length(find(rank1(1,idx.(conf.names{i})) < conf.significant));
    disp([conf.names{i} ' - sorting out: ' num2str(n_sort_out) ' sessions']);
    idx.(conf.names{i}) = idx.(conf.names{i}) & rank1(1,:) > conf.significant;
    disp([conf.names{i} ' remaining: ' num2str(length(find(idx.(conf.names{i})))) ' sessions' ]);
end

% NOTE I take pronation for sorting, because supination not for all
% available and they are strongly correlated

clear rank1 r p all_nmf all_pca all_pro all_sup h i j n_sort_out



%% Remaining Error
% remaining error is plotted only for pronation handposition because this
% is the only position which is available for all monkeys

h = figure('Visible', 'off');

% plot mean residual and mean shuffled
x       = 0:conf.max_channels;
colors  = {'b', 'g', 'c'};

for j=1:length(conf.names)
    data = vertcat(sessions(idx.(conf.names{j})).r_nmf_raw_pro);
    y    = [ones(length(find(idx.(conf.names{j}))),1)*100 data]';
    plot(x, y, 'Color', colors{j});
    hold on
end

y = [100 mean(vertcat(sessions.r_nmf_raw_pro))];
plot(x, y, 'k', 'LineWidth', 1.5);
y = [100 mean(vertcat(sessions.r_nmf_s_raw_pro))];
plot(x, y, 'r', 'LineWidth', 1.5);

% plot 15 per cent line
plot(x, ones(1,length(x)) * 15, '--k');

hold off
saveas(h, [conf.outpath  'resid_test.' conf.image_format]);
close(h);

% std of residual values for different model orders, averaged over all the
% sessions. This plot serves to see for which dimensionality reduction the
% nmf is stable
h = figure('Visible', 'off');
plot(mean(vertcat(sessions.std_nmf_raw_pro)))
title('std of resid values for different model orders');
saveas(h, [conf.outpath  'resid_test_std.' conf.image_format]);
close(h);

clear x y data map modi colors h j 



%% what is the remaining error for rank 3 model?
for i = 1:conf.n_monks
    data = vertcat(sessions(idx.(conf.names{i})).r_nmf_raw_pro);
    h = figure('Visible', 'off');
    hist(data(:,conf.dim),50);
    title(['dist of rank ' num2str(conf.dim) ' resid values -- mean ' ...
        num2str(mean(data(:, conf.dim)))]);
    saveas(h, [conf.outpath  'resid_dist_' conf.names{i} '.' conf.image_format]);
    close(h);
    disp('');
    disp(['mean of rank ' num2str(conf.dim) ' resid values: '  ...
        num2str(mean(data(:,conf.dim)))]);
end
clear data h i



%% for vega (where both handpositions are available) is there a difference
% in the residual?

for i= 1:conf.n_monks
    
    if max([sessions(idx.(conf.names{i})).hands] >1)
        h = figure('Visible', 'off');
        
        % select sessions from vega for which both handpos available
        x       = 0:conf.max_channels;
        index   = idx.(conf.names{i}) & ([sessions.hands] > 1);
        n_index = length(find(index));
        pro     = vertcat(sessions(index).r_nmf_raw_pro);
        sup     = vertcat(sessions(index).r_nmf_raw_sup);
        
        subplot 311
        plot(x, [ones(n_index,1)*100 pro]', 'b')
        title('pronation');
        
        subplot 312
        plot(x, [ones(n_index,1)*100 sup]', 'b')
        title('supination');
        
        subplot 313
        plot(x, [100 mean(pro)])
        hold on
        plot(x, [100 mean(sup)], 'k')
        hold off
        
        saveas(h, [conf.outpath  'rank1_handpos_' conf.names{i} '.' conf.image_format]);
        close(h);
        
        figure('Visible', 'off');
        subplot 221
        hist(pro(:,1), 1:100)
        title('pronation rank 1');
        
        subplot 222
        hist(pro(:,2), 1:100)
        title('pronation rank 2');
        
        subplot 223
        hist(sup(:,1), 1:100)
        title('supination rank 1');
        
        subplot 224
        hist(sup(:,2), 1:100)
        title('supination rank 2');
        
        
        saveas(h, [conf.outpath  'rank1_dist_' conf.names{i} '.' conf.image_format]);
        close(h);                
    end
end
clear x index n_index pro sup h i




%% synergy computation stuff

modi = {'_pro', '_sup'};

for i = 1:conf.n_monks

    monk = conf.names{i};
    
    if max([sessions(idx.(conf.names{i})).hands]) == 1
        res.(monk).nmf_sup         = [];
        res.(monk).synall_sup      = [];
        res.(monk).score_sup       = [];
        res.(monk).scoreall_sup    = [];
    end
        
    for j = 1:max([sessions(idx.(monk)).hands])
        
        if j == 2
            ses2take = idx.(monk) & [sessions.hands] == 2;
        else
            ses2take = idx.(monk);
        end
        
        group_nmf           = group(sessions(ses2take), ['nmf' modi{j}]);
        group_pca           = group(sessions(ses2take), ['pca' modi{j}]);
        synnmf              = group_nmf(1).center;
        synpca              = group_pca(1).center;
        
        if conf.norm_matchscore
            baseline                            = res.(monk).stats.m_base;
            [synpca, synnmf, score_group]       = match_syns(synpca, synnmf, baseline);
        else
            [synpca, synnmf, score_group]       = match_syns(synpca, synnmf);
        end
        res.(monk).(['synnmf' modi{j}])     = synnmf;
        res.(monk).(['nmfpca_sc' modi{j}])  = score_group;
        res.(monk).(['synpca' modi{j}])     = synpca;
        res.(monk).(['nmfdat'  modi{j}])    = vertcat(group_nmf.dat);
        res.(monk).(['pcadat'  modi{j}])    = vertcat(group_pca.dat);
        res.(monk).(['nmf_std' modi{j}])    = group_nmf(1).idx;
        res.(monk).(['pca_std' modi{j}])    = group_pca(1).idx;
        
        % will the synergies look the same when I compute them on the whole dat
        all = [];
        for k = find(ses2take)
            all = [all; sessions(k).mats(j).data_raw(:,res.(monk).c2take)]; %#ok<AGROW>
        end
        nmf_res = nmf_explore(all, conf);
        synall  = nmf_res.syns;

        if conf.norm_matchscore
            [~, synall, score_all]              = match_syns(synpca, synall, baseline);
        else
            [~, synall, score_all]              = match_syns(synpca, synall);
        end
        res.(monk).(['synall' modi{j}])     = synall;
        res.(monk).(['allpca_sc' modi{j}])  = score_all;
    end
end
clear synnmf synpca score_group baseline all group_nmf group_pca h i j k modi monk nmf_res score_all ses2take synall
    



%% consistency over sessions

for i = 1:conf.n_monks
    monk = conf.names{i};

    h = figure('Visible', 'off');
   
    subplot(6,4,[1 5 9]);
    imagesc((res.(monk).nmfdat_pro));
    axis off
    title(['consistency over sessions (nmf)' conf.names{i}]);
    
    subplot(6,4,[13 17 21]);
    imagesc((res.(monk).pcadat_pro));
    axis off
    title(['consistency over sessions (pca)' conf.names{i}]);
    
    subplot(6,4,2);
    imagesc(res.(monk).synall_pro);
    axis off
    title('all at once');
    
    subplot(6,4,10);
    imagesc(res.(monk).synnmf_pro);
    axis off
    title('centers (nmf)');
    
    subplot(6,4,14);
    imagesc(res.(monk).synpca_pro);
    axis off
    title('centers (pca)');
    
    p_pos = {[3 7], [11 15], [19 23]};
    for j = 1:size(res.(monk).synnmf_pro,1)
        subplot(6,4,p_pos{j})
        bar( [res.(monk).synpca_pro(j,:)' res.(monk).synnmf_pro(j,:)']);
        axis off
        text(10, 0.8, ['#' int2str(j) ' sc: ' num2str(res.(monk).nmfpca_sc_pro(j))]);
    end
    subplot(6,4,p_pos{1})
    title('pca vs. means');
    
    p_pos = {[4 8], [12 16], [20 24]};
    for j = 1:size(res.(monk).synpca_pro,1)
        subplot(6,4,p_pos{j})
        bar( [res.(monk).synpca_pro(j,:)' res.(monk).synall_pro(j,:)']);
        axis off
        text(10, 0.8, ['#' int2str(j) ' sc: ' num2str(res.(monk).allpca_sc_pro(j))]);
    end
    subplot(6,4,p_pos{1})
    title('pca vs. at once');

  
    saveas(h, [conf.outpath  'syn_consist_sessions_' conf.names{i} '.' conf.image_format]);
    close(h);
end


% stds of clustering
h = figure('Visible', 'off');

for i = 1:conf.n_monks
    subplot(conf.n_monks,2,i*2-1);
    hist(res.(conf.names{i}).nmf_std_pro);
    title([conf.names{i} ' std of clustering (nmf): ' ...
        num2str(std(hist(res.(conf.names{i}).nmf_std_pro, conf.dim)))]);
    
    subplot(conf.n_monks,2,i*2);
    hist(res.(conf.names{i}).pca_std_pro);
    title([conf.names{i} ' std of clustering (pca): ' ...
        num2str(std(hist(res.(conf.names{i}).pca_std_pro, conf.dim)))]);    
end
saveas(h, [conf.outpath  'syn_consist_sessions_std.' conf.image_format]);
close(h);


clear p_pos h i j monk




%% stability over posture (when available)

for i= 1:conf.n_monks
    
    if max([sessions(idx.(conf.names{i})).hands]) >1
        
        h = figure('Visible', 'off');    
        
        max_hands   = max([sessions(idx.(conf.names{i})).hands]);
        ses2take    = idx.(conf.names{i}) & [sessions.hands] == max_hands;
        g_pro       = group(sessions(ses2take), 'nmf_pro');
        g_sup       = group(sessions(ses2take), 'nmf_sup');
        
        pro_syns    = g_pro(1).center;
        sup_syns    = g_sup(1).center;
        
        if conf.norm_matchscore
            baseline    = res.(monk).stats.h_base;
            [pro_syns, sup_syns, scores] = match_syns(pro_syns, sup_syns, baseline);
        else
            [pro_syns, sup_syns, scores] = match_syns(pro_syns, sup_syns);
        end
        
        for j = 1:size(pro_syns,1)
            subplot(4,1,j)
            bar( [pro_syns(j,:)' sup_syns(j,:)']);
            axis off
            title(['#' int2str(j) ' sc: ' num2str(scores(j))]);
        end
        subplot(4,1,4)
        plot(pro_syns(:), sup_syns(:), '.');
        [r, p] = corrcoef(pro_syns(:), sup_syns(:));
        title(['r: ' num2str(r(2,1)) ' -- p: ' num2str(p(2,1))]);
        saveas(h, [conf.outpath  'post_consist_' conf.names{i} '.' conf.image_format]);
    end
end
clear max_hands ses2take g_pro g_sup pro_syns sup_syns p r h i j scores




        






%% synergies in PD space

for i = 1:conf.n_monks
    
    
    % roseplot 
    h = figure('Visible', 'off');
    syn      = res.(conf.names{i}).synnmf_pro;
    pds      = res.(conf.names{i}).pds(1,:);
    
    plot_rose(h, syn, pds);
    saveas(h, [conf.outpath  'syn_rose_pro_' conf.names{i} '.' conf.image_format]);
    close(h);
    
    if max([sessions(idx.(conf.names{i})).hands] >1)
        h = figure('Visible', 'off');
        syn      = res.(conf.names{i}).synnmf_sup;
        pds      = res.(conf.names{i}).pds(1,:);
        
        plot_rose(h, syn, pds);
        saveas(h, [conf.outpath  'syn_rose_sup_' conf.names{i} '.' conf.image_format]);
        close(h);
    end
end
clear h i pds syn


%% look for directedness of synergies


% different synergies
h = figure('Visible','off');
for i = 1:length(conf.names);

    % get distribution of muscle activations for comparison (bootstrap)    
    flat = res.(conf.names{i}).nmfdat_pro(:);
    pds  = res.(conf.names{i}).pds(1,:);
    
    % for all synergies
    for j = 1:conf.dim
        
        % bootstrap
        boot_dist = NaN(1, conf.n_bootstrap);
        for boot = 1:conf.n_bootstrap;
            perms           = randperm(size(flat,1));
            rand_act        = flat(perms(1:size(res.(conf.names{i}).nmfdat_pro,2)));
            boot_dist(boot) = circ_std( pds(rand_act > median(rand_act))'); 
        end
        
        
        syn     = res.(conf.names{i}).synnmf_pro(j,:);
        cstd(j) = circ_std( pds(1, syn > median(syn))'); 
    end
        
    subplot(length(conf.names), 1, i);
    [a b] = hist(boot_dist(:));
    a = a./conf.n_bootstrap;
    bar(b,a,'b');
    hold on;
    [a b] = hist(cstd);
    bar(b,a,'r');
    hold off;
    legend('shuffled data', 'synergies');
    title(conf.names{i});
end
saveas(h, [conf.outpath 'circ_std_syns.' conf.image_format]);
close(h);
clear cstd boot_dist flat pds a b boot i j perms rand_act


%% save the results
for i = 1:length(conf.names)
    nat_mov_res = res.(conf.names{i});
    save([conf.inpath 'nat_mov_res_' conf.names{i} '.mat'], 'nat_mov_res');
end
