
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
conf.names              = {'chalva', 'vega'};
conf.n_monks            = length(conf.names);
conf.significant        = 25;
conf.max_channels       = 16;
conf.dim                = 3;
conf.Niter_exploration  = 2;
conf.n_best             = 2;
conf.rank1_raw          = true;
conf.n_trials           = 50;
conf.image_format       = 'jpg';
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
        res.(conf.names{i}).stats = stats;
    else
        tmp         = load([conf.inpath 'all_data_' conf.names{i}]);
        sessions    = [sessions tmp.sessions]; %#ok<AGROW>
        res.(conf.names{i}).stats = tmp.stats;        
    end        
end

addpath('../lib'); 
addpath('../plots');
clear tmp



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

subplot(4,conf.n_monks, conf.n_monks*3+1);
hist([sessions.hands]);
title('handpositions');

saveas(h, [conf.outpath  'statistics.' conf.image_format]);
close(h);
clear y x


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
    
    disp([monk ': sort out ' num2str(length(find(out >= 4))) ...
        ' sessions because unstable pds']);
    idx.(monk)(idx.(monk)) = idx.(monk)(idx.(monk)) & out < 4;
end
clear data


%% which channels
% assign muscle mappings
     %(vega_first, vega_later, darma, chalva)
chan = {'FCU-F', 'X', 'PL-F', 'FCR-F', 'BR-B', 'X', 'FDS-F', 'FDP-F', ...
            'X', 'EDC-E', 'X', 'APL-E', 'EDC-E', 'ECR-E', 'ED45-E', 'ECU-E';
        'BR-B', 'EDC-E', 'APL-E', 'ECU-E', 'FCR-F', 'APL-E', 'ED45-E', 'ED23-E', ...
            'ECU-E', 'BR-B', 'PL-F', 'FCR-F', 'X', 'X', 'X', 'FDS-F';
        'ECU-E', 'ED45-E', 'EDC-E', 'APL-E', 'ECR-E', 'ED23-E', 'BIC-P', 'BIC-P', ...
            'FDS-F', 'PL-F', 'FCU-F', 'FCR-F', 'PT-F', 'FDP-F', 'TRIC-P', 'BIC-P';
        'FCU-F', 'FDS-F', 'PL-F', 'FCR-F', 'PT-F', 'FDP-F', 'PL-F', 'BIC-P', ...
            'ECU-E', 'EDC-E', 'ED45-E', 'ECR-E', 'ED23-E', 'APL-E', 'ECR-E', 'TRIC-P'};
        
tmp = intersect({chan{1,:}},{chan{2,:}});
tmp = intersect(tmp,{chan{3,:}});
disp('');
disp([num2str(length(tmp)) ' channels in common for all monkeys']);

% only use channels which are available for all sessions of a monk
for i = 1:conf.n_monks
    
    all_chan                     = vertcat(sessions(idx.(conf.names{i})).channels);
    res.(conf.names{i}).c2take   = all(all_chan);
end

clear tmp chan all_chan




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

clear rank1 r p all_nmf all_pca all_pro all_sup





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

clear x y data map modi


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
clear data






%% for vega (where both handpositions are available) is there a difference
% in the residual?

for i= 1:conf.n_monks
    
    if max([sessions(idx.(conf.names{i})).hands] >1)
        figure('Visible', 'off');
        
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
                
        clear x index n_index pro sup 
    end
end



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
        res.(monk).nmfdat   = vertcat(group_nmf.dat);
        res.(monk).pcadat   = vertcat(group_pca.dat);
        res.(monk).nmf_std  = group_nmf(1).idx;
        res.(monk).pca_std  = group_pca(1).idx;
        
        baseline                            = res.(monk).stats.m_base;
        [synpca, synnmf, score_group]       = match_syns(synpca, synnmf, baseline);
        res.(monk).(['synnmf' modi{j}])     = synnmf;
        res.(monk).(['nmfpca_sc' modi{j}])  = score_group;
        res.(monk).(['synpca' modi{j}])     = synpca;
        
        % will the synergies look the same when I compute them on the whole dat
        all = [];
        for k = find(ses2take)
            all = [all; sessions(k).mats(j).data_raw(:,res.(monk).c2take)]; %#ok<AGROW>
        end
        nmf_res = nmf_explore(all, conf);
        synall  = nmf_res.syns;
        
        [~, synall, score_all]              = match_syns(synpca, synall, baseline);
        res.(monk).(['synall' modi{j}])     = synall;
        res.(monk).(['allpca_sc' modi{j}])  = score_all;
    end
end
clear synnmf synpca score_group all baseline
    



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
    hist(res.(conf.names{i}).nmf_std);
    title([conf.names{i} ' std of clustering (nmf): ' ...
        num2str(std(hist(res.(conf.names{i}).nmf_std, conf.dim)))]);
    
    subplot(conf.n_monks,2,i*2);
    hist(res.(conf.names{i}).pca_std);
    title([conf.names{i} ' std of clustering (pca): ' ...
        num2str(std(hist(res.(conf.names{i}).pca_std, conf.dim)))]);    
end
saveas(h, [conf.outpath  'syn_consist_sessions_std.' conf.image_format]);
close(h);


clear p_pos







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
        
        baseline    = res.(monk).stats.h_base;
        [pro_syns, sup_syns, scores] = match_syns(pro_syns, sup_syns, baseline);
        
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
clear max_hands ses2take g_pro g_sup pro_syns sup_syns p r







%% stability of prefered directions (over sessions)
for i = 1:conf.n_monks

    max_hands   = max([sessions(idx.(conf.names{i})).hands]);
    ses2take    = idx.(conf.names{i}) & [sessions.hands] == max_hands;
    n2take      = length(find(ses2take));
    c2take_idx  = find(res.(conf.names{i}).c2take);
    
    
    h = figure('Visible', 'off');    
    all = vertcat(sessions(ses2take).pd);            
    for k = 1:length(c2take_idx)
        
        subplot(ceil(length(c2take_idx)/2),2,k);
        if max_hands == 1
            [x y] = pol2cart(all(:,c2take_idx(k)), ones(size(all(:,c2take_idx(k)))));
            feather(x, y, 'b');
            [~, res.(conf.names{i}).cstd]  = circ_std(all);
            res.(conf.names{i}).pds        = circ_mean(all);
        else
            id_pro = 1:2:2*n2take-1;
            id_sup = 2:2:2*n2take;
            [x y] = pol2cart(all(id_pro,c2take_idx(k)), ones(size(all(id_pro,c2take_idx(k)))));
            feather(x, y, 'b');
            hold on
            [x y] = pol2cart(all(id_sup,c2take_idx(k)), ones(size(all(id_sup,c2take_idx(k)))));
            feather(x, y, 'r');
            hold off
            [~, res.(conf.names{i}).cstd(1,:)] = circ_std(all(id_pro,:));
            [~, res.(conf.names{i}).cstd(2,:)] = circ_std(all(id_sup,:));
            res.(conf.names{i}).pds(1,:)       = circ_mean(all(id_pro,:));
            res.(conf.names{i}).pds(2,:)       = circ_mean(all(id_sup,:));
        end            
    end
    saveas(h, [conf.outpath  'pd_consist_feather_' conf.names{i} '.' conf.image_format]);
    close(h);
    
    
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
                res.(conf.names{i}).cstd(j,c2take_idx), 'UniformOutput', false))
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


clear x y all c2take_idx n_hands monk_first col in_deg inds id_pro id_sup
        






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
    
    nat_mov_res = res.(conf.names{i});
    save([conf.inpath 'nat_mov_res_' conf.names{i} '.mat'], 'nat_mov_res');
    
end



clear h_fake pds rose_agg syn



