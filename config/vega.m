

config = struct;
config.create_roseplots = 0;
config.nmf_stab_plots   = 0;
config.fetch_data       = 1;
config.Niter_exploration= 50;
config.n_iter_stab      = 50;
config.Niter_res_test   = 20;
config.dimensions       = 3;
config.significant      = 25;
config.channels2take    = find(logical([0  1  1  1  1  1  1  1  1  1  1  1  0  0  0  0 ]));
config.n_used_channels  = length(config.channels2take);
config.emgdat_path      = '/Volumes/dedan/yifat_lab_data/vega/EMGdat/';
config.outpath          = '~/Documents/uni/yifat_lab/results/nonevoked_syns/';
config.image_format     = 'pdf';
config.res_folder       = '~/Documents/uni/yifat_lab/results/';
config.pd_folder        = '/Volumes/DEDAN_DISK/vega/pd_files/';
config.names            = {'_pro', '_sup',''};
config.opt              = statset('MaxIter',5);