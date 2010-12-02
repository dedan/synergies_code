
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
conf.inpath             = '~/Documents/uni/yifat_lab/results/data/';   
conf.names              = {'chalva'};
conf.n_monks            = length(conf.names);
conf.significant        = 25;
conf.max_channels       = 16;
conf.dimensions         = 3;
conf.Niter_exploration  = 5;
conf.n_best             = 5;
conf.rank1_raw          = true;
conf.image_format       = 'jpg';
mymap = [linspace(145/255,178/255,32)' linspace(167/255,213/255,32)' linspace(216/255,111/255,32)'];
mymap = [mymap; linspace(178/255,251/255,32)' linspace(213/255,147/255,32)' linspace(111/255,24/255,32)'];
conf.map = mymap;

diary([conf.outpath 'log.txt']);

if conf.n_best > conf.Niter_exploration
    disp('n_best has to be smaller than Niter_exploration');
end
clear mymap


%% load data
load([conf.inpath 'all_data_dummy']);
addpath('../lib'); 

% sort out the first vega sessions because they were recorded from
% different muscles and have only 8 muscles in common with the later
% sessions which is I think not enough
not_first_vega     = ~strcmp('vega', {sessions.monk}) | [sessions.id] > 26;
for i = find(~not_first_vega)
    sessions(i).monk = 'first_vega'; %#ok<SAGROW>
end
clear not_first_vega


%% statistics
% First lets get some statistics on the Data we have available

for i = 1:length(sessions)
    sessions(i).target1     = size(sessions(i).mats(1).data,1); %#ok<SAGROW>
    sessions(i).n_channels  = length(find(sessions(i).channels)); %#ok<SAGROW>
    
    if length(sessions(i).mats) == 2
        sessions(i).target2    = size(sessions(i).mats(2).data,1); %#ok<SAGROW>
    else
        sessions(i).target2    = -1; %#ok<SAGROW>
    end
end

h = figure('Visible', 'off');
disp(['data from in total: ' num2str(length(sessions)) ' sessions']);

for i = 1:conf.n_monks
    idx.(conf.names{i})    = strcmp(conf.names{i}, {sessions.monk});
    
    disp([conf.names{i} ': ' num2str(length(find(idx.(conf.names{i})))) ' sessions' ]);

    % distribution of targets
    subplot(3,conf.n_monks,i);
    
    [y x] = hist([sessions(idx.(conf.names{i})).target1],-1:8);
    bar(x,y,'r');
    
    hold on
    [y x] = hist([sessions(idx.(conf.names{i})).target2],-1:8);
    bar(x+0.2,y,'b');
    hold off
    title(conf.names{i});
    xlabel('targets');

    
    % distribution of used channels
    subplot(3, conf.n_monks, conf.n_monks +i);
    
    [y x] = hist([sessions(idx.(conf.names{i})).n_channels],1:16);
    bar(x,y,'b');
    xlabel('channels');
end

subplot(3,conf.n_monks, conf.n_monks*2+1);
hist([sessions.hands]);
title('handpositions');

saveas(h, [conf.outpath  'statistics.' conf.image_format]);
close(h);
clear y x




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

% which channels are available for all sessions of a monkey
for i = 1:conf.n_monks
    
    % only use channels which are available for all sessions of a monk
    all_chan = vertcat(sessions(idx.(conf.names{i})).channels);
    c2take   = all(all_chan);
    for j = find(idx.(conf.names{i}))
        sessions(j).c2take = c2take; %#ok<SAGROW>
    end
end

clear tmp chan all_chan c2take




%% rank1 analysis

% plot the rank1 values for all sessions and the different tests to see
% relation between the different tests


% first store them in a convenient matrix
rank1 = NaN(2,length(sessions));
for i = 1:length(sessions)
    if conf.rank1_raw
        rank1(1,i)  = sessions(i).r_nmf_raw_pro(1);
        rank1(2,i)  = sessions(i).r_pca_raw_pro(1);
    else
        rank1(1,i)  = sessions(i).r_nmf_pro(1);
        rank1(2,i)  = sessions(i).r_pca_pro(1);
    end
end

h = figure('Visible', 'off');

% plot the rank1 values for pca and nmf and also the difference between
% them. furthermore the line at the value at which sessions are chosen as
% significant
plot(rank1','.');
hold on;
plot(ones(1,length(rank1)) * conf.significant);

% seperating lines between the monkeys
for i=2:conf.n_monks-1
    plot(ones(1,10)*find(idx.(conf.names{i}), 1, 'last' ), 1:10:100, 'r*');
end
hold off;
legend('nmf', 'pca', 'Location', 'NorthWest');
for i = 1:conf.n_monks
    text(find(idx.(conf.names{i}), 1 )+2, 70, conf.names{i});
end

% also calculate the correlation of the two rank1 values
[r p] = corrcoef(rank1');
title(['nmf vs. pcaica, r: ' num2str(r(1,2)) ' - p: ' num2str(p(1,2))]);

if conf.rank1_raw
    saveas(h, [conf.outpath  'rank1_raw.' conf.image_format]);
else
    saveas(h, [conf.outpath  'rank1.' conf.image_format]);
end
    
close(h);


% sort out boring sessions
% we sort out sessions where the rank1 model already has an remaining error
% lower then 25 %
disp('')
disp('rank1 filtering');
disp(['sorting out: ' num2str(length(find(rank1(1,:) < conf.significant))) ' sessions']);
for i = 1:conf.n_monks
    idx.(conf.names{i}) = idx.(conf.names{i}) & rank1(1,:) > conf.significant;
    disp([conf.names{i} ' remaining: ' num2str(length(find(idx.(conf.names{i})))) ' sessions' ]);
end



clear rank1 r p








%% Remaining Error
% remaining error is plotted only for pronation handposition because this
% is the only position which is available for all monkeys

% make all test fields the same length (easier access later on)


h = figure('Visible', 'off');

% plot mean residual and mean shuffled
x = 0:conf.max_channels;
modi = {'_raw', ''};
map = jet;

for i = 1:2
    subplot(2,1,i);
    for j=1:length(conf.names)
        data = vertcat(sessions(idx.(conf.names{j})).(['r_nmf' modi{i} '_pro']));
        y    = [ones(length(find(idx.(conf.names{j}))),1)*100 data]';
        plot(x, y, 'Color', map(j*15,:));
        hold on
    end
    y = [100 mean(vertcat(sessions.(['r_nmf' modi{i} '_pro'])))];
    plot(x, y, 'k', 'LineWidth', 1.5);
    y = [100 mean(vertcat(sessions.(['r_nmf_s' modi{i} '_pro'])))];
    plot(x, y, 'r', 'LineWidth', 1.5);
    hold off
end
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




%% for vega (where both handpositions are available) is there a difference
% in the residual?

% TODO nicht nur fuer vega, auch fuer darma haben wir zwei
% TODO und auch fuer first_vega

if any(strcmp(conf.names, 'vega'))
    figure('Visible', 'off');
    
    % select sessions from vega for which both handpos available
    x       = 0:conf.max_channels;
    index   = idx.vega & ([sessions.hands] > 1);
    n_index = length(find(index));
    pro     = vertcat(sessions(index).r_nmf_raw_pro);
    sup     = vertcat(sessions(index).r_nmf_raw_sup);
    
    subplot 311
    plot(x, [ones(n_index,1)*100 pro]', 'b')
    
    subplot 312
    plot(x, [ones(n_index,1)*100 sup]', 'b')
    
    subplot 313
    plot(x, [100 mean(pro)])
    hold on
    plot(x, [100 mean(sup)], 'k')
    hold off
    
    saveas(h, [conf.outpath  'rank1_handpos.' conf.image_format]);
    close(h);
    
    figure('Visible', 'off');
    subplot 221
    hist(pro(:,1))
    title('pronation rank 1');
    
    subplot 222
    hist(pro(:,2))
    title('pronation rank 2');
    
    subplot 223
    hist(sup(:,1))
    title('supination rank 1');
    
    subplot 224
    hist(sup(:,2))
    title('supination rank 2');
    
    
    saveas(h, [conf.outpath  'rank12_dist.' conf.image_format]);
    close(h);
    
    [h, p] = kstest2(pro(:,1), sup(:,1));
    disp('similarity of rank 1 distributions for pronation and supination');
    if h == 1
        disp(['rank 1 distributions not similar, p: ' num2str(p)]);
    else
        disp(['rank 1 distributions are similar, p: ' num2str(p)]);
    end
    
    clear x index n_index pro sup h p
end




%% compute the synergies 
% all is done only by using nmf. that it leads to the same results as
% pcaica will be shown elsewhere

if exist([conf.inpath 'all_data_syn.mat'], 'file')
    load([conf.inpath 'all_data_syn']);
else
    
    for i = 1:conf.n_monks
        
        for j = find(idx.(conf.names{i}))
            
            disp(['monk: ' conf.names{i} ' session: ' num2str(j)]);
            data    = sessions(j).mats(1).data_raw(:,sessions(j).c2take);
            nmf_res = nmf_explore(data, conf);
            sessions(j).nmf_pro     = nmf_res.syns; %#ok<SAGROW>
            sessions(j).nmf_pro_std = nmf_res.std; %#ok<SAGROW>
            sessions(j).pca_pro     = pcaica(data, conf.dimensions)'; %#ok<SAGROW>
            
            
            if length(sessions(j).mats) > 1
                data    = sessions(j).mats(2).data_raw(:,sessions(j).c2take);
                nmf_res = nmf_explore(data, conf);
                sessions(i).nmf_sup     = nmf_res.syns; %#ok<SAGROW>
                sessions(i).nmf_sup_std = nmf_res.std; %#ok<SAGROW>
                sessions(j).pca_sup     = pcaica(data, conf.dimensions)'; %#ok<SAGROW>
            end
        end
    end
    save([conf.inpath 'all_data_syn'], 'sessions');
end

clear all_chan data nmf_res data




%% synergy analysis



% consistency over sessions
for i = 1:conf.n_monks

    h = figure('Visible', 'off');
    
    group_nmf = group(sessions(idx.(conf.names{i})), 'nmf_pro');    
    subplot(6,4,[1 5 9]);
    imagesc(vertcat(group_nmf.dat));
    axis off
    title(['consistency over sessions (nmf)' conf.names{i}]);
    
    group_pca = group(sessions(idx.(conf.names{i})), 'pca_pro');    
    subplot(6,4,[13 17 21]);
    imagesc(vertcat(group_pca.dat));
    axis off
    title(['consistency over sessions (pca)' conf.names{i}]);
    
    % will the synergies look the same when I compute them on the whole dat
    all = [];
    for k = find(idx.(conf.names{i}))
        all = [all; sessions(k).mats(1).data_raw(:,sessions(k).c2take)]; %#ok<AGROW>
    end
    nmf_res = nmf_explore(all, conf);
    res     = nmf_res.syns;
    subplot(6,4,2);
    imagesc(nmf_res.syns);
    axis off
    title('all at once');
    
    subplot(6,4,10);
    imagesc(group_nmf(1).center);
    axis off
    title('centers (nmf)');
    stds_n{i} = group_nmf(1).idx;   %#ok<SAGROW>
    
    subplot(6,4,14);
    imagesc(group_pca(1).center);
    axis off
    title('centers (pca)');
    stds_p{i} = group_pca(1).idx;   %#ok<SAGROW>

    
    synnmf = normr(group_nmf(1).center);
    synpca = normr(group_pca(1).center);
    synall = normr(res);
    
    [scores ind] = matchNscore(synpca', synnmf');
    synnmf       = synnmf(ind,:);
    
    p_pos = {[3 7], [11 15], [19 23]};
    for j = 1:size(synnmf,1)
        subplot(6,4,p_pos{j})
        bar( [synpca(j,:)' synnmf(j,:)']);
        axis off
        title(['#' int2str(j) ' sc: ' num2str(scores(j))]);
    end
    [r p] = corrcoef([synpca(:) synnmf(:)]);
    disp(['nmf vs. pca, r: ' num2str(r(1,2)) ' - p: ' num2str(p(1,2))]);

    
    [scores ind] = matchNscore(synpca', synall');
    synall       = synall(ind,:);
    
    p_pos = {[4 8], [12 16], [20 24]};
    for j = 1:size(synpca,1)
        subplot(6,4,p_pos{j})
        bar( [synpca(j,:)' synall(j,:)']);
        axis off
        title(['#' int2str(j) ' sc: ' num2str(scores(j))]);
    end
    [r p] = corrcoef([synpca(:) synall(:)]);
    disp(['nmf vs. pca, r: ' num2str(r(1,2)) ' - p: ' num2str(p(1,2))]);
  
    saveas(h, [conf.outpath  'syn_consist_sessions_' conf.names{i} '.' conf.image_format]);
    close(h);
end


% stds of clustering
h = figure('Visible', 'off');

for i = 1:conf.n_monks
    subplot(conf.n_monks,2,i*2-1);
    hist(stds_n{i});
    title([conf.names{i} ' std of clustering (nmf): ' num2str(std(hist(stds_n{i}, conf.dimensions)))]);
    
    subplot(conf.n_monks,2,i*2);
    hist(stds_p{i});
    title([conf.names{i} ' std of clustering (pca): ' num2str(std(hist(stds_p{i}, conf.dimensions)))]);    
end
saveas(h, [conf.outpath  'syn_consist_sessions_std.' conf.image_format]);
close(h);


clear flat grouped stds_n stds_p nmf_res all synnmf synpca synall






% stability over posture (when available)
% stability over monkeys 



%% stability of prefered directions (over sessions)

for i = 1:conf.n_monks

    n_hands = sessions(idx.(conf.names{i})(1)).hands;
    c2take  = find(sessions(idx.(conf.names{i})(1)).c2take);
    
    h = figure('Visible', 'off');
    for j = 1:n_hands
        subplot(n_hands,1,j)
        for k = find(idx.(conf.names{i}))
            col = [k k k]/length(idx.(conf.names{i}));
            plot(c2take,sessions(k).pd_deg(j,c2take),'^','MarkerSize',10, 'MarkerFaceColor', col);
            hold on
        end
        grid
        hold off
    end
    saveas(h, [conf.outpath  'pd_consist_' conf.names{i} '.' conf.image_format]);
    close(h);
    
    
    h = figure('Visible', 'off');
        
    all = vertcat(sessions(idx.(conf.names{i})).pd);
            
    for k = c2take
        subplot(length(c2take)/2,2,k);
        if n_hands == 1
            [x y] = pol2cart(all(:,k), ones(size(all(:,k))));
            feather(x, y, 'b');
        else
            [x y] = pol2cart(all(1:2:length(c2take)-1,k), ones(size(all(1:2:length(c2take)-1,k))));
            feather(x, y, 'b');
            hold on
            [x y] = pol2cart(all(2:2:length(c2take),k), ones(size(all(2:2:length(c2take),k))));
            feather(x, y, 'r');
            hold off
            legend('pronation', 'supination');
        end            
    end
    saveas(h, [conf.outpath  'pd_consist_feather' conf.names{i} '.' conf.image_format]);
    close(h);
end

clear x y all c2take n_hands
        


















