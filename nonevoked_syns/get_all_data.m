
% This is a runscript to do all the computationally expensive stuff like
% residual tests. Once this is run, the files which are created can be used
% for further analysis

clear

addpath('../lib'); %#ok<MCAP>

% general settings
config = struct;
config.Niter_exploration= 3;
config.n_iter_stab      = 3;
config.Niter_res_test   = 1;
config.opt              = statset('MaxIter',5);
config.outpath          = '~/Documents/uni/yifat_lab/results/natural_mov/';
config.res_folder       = '~/Documents/uni/yifat_lab/results/';

diary([config.outpath 'log.txt']);


% for chalva
config.monk             = 'chalva';
config.emgdat_path      = ['/Volumes/LAB/' config.monk filesep 'EMGdat' filesep];
config.pd_folder        = ['/Volumes/LAB/' config.monk filesep 'pd_files' filesep];
config.names            = {'_pro'};

sessions = nonevoked_sessions(config);



% for vega
config.monk             = 'vega';
config.emgdat_path      = ['/Volumes/LAB/' config.monk filesep 'EMGdat' filesep];
config.pd_folder        = ['/Volumes/LAB/' config.monk filesep 'pd_files' filesep];
config.names            = {'_pro', '_sup', ''};

sessions = [sessions nonevoked_sessions(config)];


% for darma
config.monk             = 'darma';
config.emgdat_path      = ['/Volumes/LAB/' config.monk filesep 'EMGdat' filesep];
config.pd_folder        = ['/Volumes/LAB/' config.monk filesep 'pd_files' filesep];
config.names            = {'_pro', '_sup', ''};

sessions = [sessions nonevoked_sessions(config)];


% assign muscle mappings
vega_first = {'FCU-F', 'X', 'PL-F', 'FCR-F', 'BR-B', 'X', 'FDS-F', 'FDP-F', ...
 'X', 'EDC-E', 'X', 'APL-E', 'EDC-E', 'ECR-E', 'ED45-E', 'ECU-E'};
vega_later = {'BR-B', 'EDC-E', 'APL-E', 'ECU-E', 'FCR-F', 'APL-E', 'ED45-E', 'ED23-E', ...
    'ECU-E', 'BR-B', 'PL-F', 'FCR-F', 'X', 'X', 'X', 'FDS-F'};
darma = {'ECU-E', 'ED45-E', 'EDC-E', 'APL-E', 'ECR-E', 'ED23-E', 'BIC-P', 'BIC-P', ...
    'FDS-F', 'PL-F', 'FCU-F', 'FCR-F', 'PT-F', 'FDP-F', 'TRIC-P', 'BIC-P'};
chalva = {'FCU-F', 'FDS-F', 'PL-F', 'FCR-F', 'PT-F', 'FDP-F', 'PL-F', 'BIC-P', ...
    'ECU-E', 'EDC-E', 'ED45-E', 'ECR-E', 'ED23-E', 'APL-E', 'ECR-E', 'TRIC-P'};


save([config.outpath 'all_data'], 'sessions');
diary('off');
disp('finished');
