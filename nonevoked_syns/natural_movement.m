
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
conf.res_folder         = '~/Documents/uni/yifat_lab/results/';   
conf.names              = {'vega', 'chalva', 'darma'};
conf.n_monks            = length(conf.names);
conf.significant        = 25;
conf.max_channels       = 16;
conf.dimensions         = 3;
conf.Niter_exploration  = 5;
conf.n_best             = 5;
mymap = [linspace(145/255,178/255,32)' linspace(167/255,213/255,32)' linspace(216/255,111/255,32)'];
mymap = [mymap; linspace(178/255,251/255,32)' linspace(213/255,147/255,32)' linspace(111/255,24/255,32)'];
conf.map = mymap;



diary([conf.outpath 'log.txt']);

if conf.n_best > conf.Niter_exploration
    disp('das geht doch nicht');
end


%% load data
load([conf.outpath 'all_data_dummy']);
addpath('../lib'); 


%% statistics
% First lets get some statistics on the Data we have available

stats = struct;
for i = 1:length(sessions)
    stats(i).monk           = sessions(i).monk;
    stats(i).target1        = size(sessions(i).mats(1).data,1);
    stats(i).channels       = size(sessions(i).mats(1).data,2);
    stats(i).hands          = sessions(i).hands;
    
    if length(sessions(i).mats) == 2
        stats(i).target2    = size(sessions(i).mats(2).data,1);
    else
        stats(i).target2    = -1;
    end
end

figure(1);
disp(['data from in total: ' num2str(length(sessions)) ' sessions']);

for i = 1:conf.n_monks
    idx.(conf.names{i})    = strmatch(conf.names{i}, {sessions.monk});
    
    disp([conf.names{i} ': ' num2str(length(idx.(conf.names{i}))) ' sessions' ]);

    % distribution of targets
    subplot(3,conf.n_monks,i);
    
    [y x] = hist([stats(idx.(conf.names{i})).target1],-1:8);
    bar(x,y,'r');
    
    hold on
    [y x] = hist([stats(idx.(conf.names{i})).target2],-1:8);
    bar(x+0.2,y,'b');
    hold off
    title(conf.names{i});
    xlabel('targets');

    
    % distribution of used channels
    subplot(3, conf.n_monks, conf.n_monks +i);
    
    [y x] = hist([stats(idx.(conf.names{i})).channels],1:16);
    bar(x,y,'b');
    xlabel('channels');
end

subplot(3,conf.n_monks, conf.n_monks*2+1);
hist([sessions.hands]);
title('handpositions');



%% rank1 analysis

% plot the rank1 values for all sessions and the different tests to see
% relation between the different tests


% first store them in a convenient matrix
rank1 = NaN(2,length(sessions));
for i = 1:length(sessions)
    rank1(1,i)  = sessions(i).r_nmf_pro(1);
    rank1(2,i)  = sessions(i).r_pca_pro(1);
end

h = figure(2);

% plot the rank1 values for pca and nmf and also the difference between
% them. furthermore the line at the value at which sessions are chosen as
% significant
plot(rank1','.');
hold on;
plot(ones(1,length(rank1)) * conf.significant);

% seperating lines between the monkeys
for i=1:conf.n_monks-1
    plot(ones(1,10)*max(idx.(conf.names{i})), 1:10:100, 'r*');
end
hold off;
legend('nmf', 'pca', 'Location', 'NorthWest');
text(min(idx.(conf.names{1}))+2, 70, conf.names{1});
text(min(idx.(conf.names{2}))+2, 70, conf.names{2});
text(min(idx.(conf.names{3}))+2, 70, conf.names{3});


% also calculate the correlation of the two rank1 values
[r p] = corrcoef(rank1');
title(['nmf vs. pcaica, r: ' num2str(r(1,2)) ' - p: ' num2str(p(1,2))]);

%saveas(h, [config.outpath  'rank1.' config.image_format]);
%close(h);
clear rank1





%% Remaining Error
% remaining error is plotted only for pronation handposition because this
% is the only position which is available for all monkeys

% make all test fields the same length (easier access later on)


h = figure(3);

% plot mean residual and mean shuffled
x = 0:conf.max_channels;
modi = {'_raw', ''};
map = jet;

for i = 1:2
    subplot(2,1,i);
    for j=1:length(conf.names)
        data = vertcat(sessions(idx.(conf.names{j})).(['r_nmf' modi{i} '_pro']));
        y    = [ones(length(idx.(conf.names{j})),1)*100 data]';
        plot(x, y, 'Color', map(j*20,:));
        hold on
    end
    y = [100 mean(vertcat(sessions.(['r_nmf' modi{i} '_pro'])))];
    plot(x, y, 'r', 'LineWidth', 1.5);
    y = [100 mean(vertcat(sessions.(['r_nmf_s' modi{i} '_pro'])))];
    plot(x, y, 'k', 'LineWidth', 1.5);
    hold off
end


%clear x y data map modi




%% for vega (where both handpositions are available) is there a difference
% in the residual?

% TODO nicht nur fuer vega, auch fuer darma haben wir zwei

figure(4)

% select sessions from vega for which both handpos available
x       = 0:conf.max_channels;
index   = [sessions(idx.vega).hands] > 1;
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

% TODO what is the correct test for this
%ranksum(pro(:,1), sup(:,1))




%% compute the synergies 
% all is done only by using nmf. that it leads to the same results as
% pcaica will be shown elsewhere

if exist([conf.outpath 'all_data_syn.mat'], 'file')
    load([conf.outpath 'all_data_syn']);
else
    
    for i = 1:length(sessions)
        
        disp(['computing synergies for session: ' num2str(i)]);
        nmf_res = nmf_explore(sessions(i).mats(1).data_raw, conf);
        sessions(i).syn_pro = nmf_res.syns; %#ok<SAGROW>
        
        if length(sessions(i).mats) > 1
            nmf_res = nmf_explore(sessions(i).mats(2).data_raw, conf);
            sessions(i).syn_sup = nmf_res.syns; %#ok<SAGROW>
        end
    end
    save([conf.outpath 'all_data_syn'], 'sessions');
end



%% synergy analysis



% consistency over sessions

% nmf seems to be quite stable, nevertheless I take the center of the
% prototypes found in the five runs

for i = 1:conf.n_monks
    
    grouped = group(sessions(idx.(conf.names{i})), 'syn_pro');
    fin_syns.(['nmf' config.names{i}]) = grouped.center;
    
    flat = vertcat(grouped.dat);
        
    h = figure(5);
    subplot(4,1,1:3);
    imagesc(flat);
    colormap(conf.map);
    axis off
    title(['consistency of synergists over sessions ' conf.names{i}]);
    
    subplot(4,1,4);
    imagesc(grouped(1).center);
    colormap(conf.map);
    axis off
    title(['centers ' config.names{i}]);
%     saveas(h, [config.outpath  'syn_consist_sessions_' config.names{i} '.' config.image_format]);
%     close(h);
    
    stds(i,:) = grouped(1).idx;   %#ok<SAGROW>
end

h = figure(6);
for i = 1:length(conf.names)
    subplot(3,1,i);
    hist(stds(i,:));
    title([conf.names{i} ' std of clustering: ' num2str(std(hist(stds(i,:), conf.dimensions)))]);
end
% saveas(h, [config.outpath  'syn_consist_sessions_std.' config.image_format]);
% close(h);


% clear flat grouped all_names






% stability over posture (when available)
% stability over monkeys 


        


















