
% This is a runscript to do all the computationally expensive stuff like
% residual tests. Once this is run, the files which are created can be used
% for further analysis

clear

addpath('../lib'); %#ok<MCAP>

% general settings
config = struct;
config.Niter_exploration= 2;
config.n_iter_stab      = 2;
config.Niter_res_test   = 2;
config.opt              = statset('MaxIter',5);
config.max_channels     = 16;
config.outpath          = 'E:\results\natural_mov\';
config.res_folder       = 'E:\results\';
% config.outpath          = '~/Documents/uni/yifat_lab/results/natural_mov/';
% config.res_folder       = '~/Documents/uni/yifat_lab/results/';
config.data_path        = 'E:\';

diary([config.outpath 'log.txt']);

% for vega
config.monk             = 'vega';
config.emgdat_path      = [config.data_path config.monk filesep 'EMGdat' filesep];
config.pd_folder        = [config.data_path config.monk filesep 'pd_files' filesep];
config.names            = {'_pro', '_sup', ''};

sessions = nonevoked_sessions(config);

% for chalva
config.monk             = 'chalva';
config.emgdat_path      = [config.data_path config.monk filesep 'EMGdat' filesep];
config.pd_folder        = [config.data_path config.monk filesep 'pd_files' filesep];
config.names            = {'_pro'};

sessions = [sessions nonevoked_sessions(config)];


% for darma
config.monk             = 'darma';
config.emgdat_path      = [config.data_path config.monk filesep 'EMGdat' filesep];
config.pd_folder        = [config.data_path config.monk filesep 'pd_files' filesep];
config.names            = {'_pro', '_sup', ''};

sessions = [sessions nonevoked_sessions(config)];





save([config.outpath 'all_data'], 'sessions');
diary('off');
disp('finished');
