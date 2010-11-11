

config = struct;
config.create_roseplots = 0;
config.nmf_stab_plots   = 0;
config.fetch_data       = 1;
config.Niter_exploration= 50;
config.n_iter_stab      = 50;
config.Niter_res_test   = 5;
config.dimensions       = 3;
config.significant      = 25;
config.normalization    = 1;
% normalization 1: mean activity in window after TO - mean activity before
% normalization 2: mean activity in window after TO / mean activity before
config.channels2take    = find(logical([1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1 ]));
config.n_used_channels  = length(config.channels2take);
config.monk             = 'chalva';
config.emgdat_path      = ['/Volumes/LAB/' config.monk filesep 'EMGdat' filesep];
config.outpath          = '~/Documents/uni/yifat_lab/results/nonevoked_syns/';
config.image_format     = 'pdf';
config.res_folder       = '~/Documents/uni/yifat_lab/results/';
config.pd_folder        = ['/Volumes/LAB/' config.monk filesep 'pd_files' filesep];
config.names            = {'_pro', '_sup',''};
config.opt              = statset('MaxIter',5);