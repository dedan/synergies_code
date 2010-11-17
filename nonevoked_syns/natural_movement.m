
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

%% confurations

clear;
clc;
conf = struct;
conf.opt              = statset('MaxIter',5);
conf.outpath          = '~/Documents/uni/yifat_lab/results/natural_mov/';
conf.res_folder       = '~/Documents/uni/yifat_lab/results/';   
conf.names            = {'vega', 'darma', 'chalva'};
conf.n_monks          = length(conf.names);

diary([conf.outpath 'log.txt']);


%% load data
load([conf.outpath 'all_data_dummy']);


%% statistics
% First lets get some statistics on the Data we have available

stats = struct;
for i = 1:length(sessions)
    stats(i).monk           = sessions(i).monk;
    stats(i).target1        = size(sessions(i).mats(1).data,1);
    stats(i).channels       = size(sessions(i).mats(1).data,2);
    stats(i).hands          = 1;
    
    if length(sessions(i).mats) == 3
        stats(i).target2    = size(sessions(i).mats(2).data,1);
        stats(i).hands      = 2;
    else
        stats(i).target2    = -1;
    end
end


figure(1);
for i = 1:conf.n_monks
    idx.(conf.names{i})    = strmatch(conf.names{i}, {sessions.monk});

    % distribution of targets
    subplot(2,conf.n_monks,i);
    
    [y x] = hist([stats(idx.(conf.names{i})).target1],-1:8);
    bar(x,y,'r');
    
    hold on
    [y x] = hist([stats(idx.(conf.names{i})).target2],-1:8);
    bar(x+0.2,y,'b');
    hold off
    title(conf.names{i});
    xlabel('targets');

    
    % distribution of used channels
    subplot(2, conf.n_monks, conf.n_monks +i);
    
    [y x] = hist([stats(idx.(conf.names{i})).channels],1:16);
    bar(x,y,'b');
    xlabel('channels');
end




%% Remaining Error




