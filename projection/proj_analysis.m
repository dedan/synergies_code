
% analysis similarity between results, originally thought to be used only
% to see wether there is a similarity between the evoked responses and the
% synergies found during natural movement. But can also be used, for
% example to look wether pronation or supination data has stronger
% influence on the synergies, found from combined data (pro-sup together)


conf.monks          = {'chalva', 'vega'};
conf.res_folder     = '/Volumes/LAB/results/';
conf.out            = '~/Documents/uni/yifat_lab/results/';
conf.n_boot         = 10000;
conf.noise          = 0.5;
conf.image_format   = 'pdf';

resps       = struct([]);


% load a result file from the syn_analysis
for i = 1:length(conf.monks)
    monk = conf.monks{i};

    add_mapinfo([conf.res_folder 'data' filesep], {monk})

    % load the nonevoked results
    tmp     = load([conf.res_folder 'data' filesep 'nat_mov_res_' monk '.mat']);
    nat_mov_res.(monk) = tmp.nat_mov_res;

    if isempty(resps)
        load([conf.res_folder 'data' filesep 'evoked_data_map_' monk '.mat']);
    else
        tmp     = load([conf.res_folder 'data' filesep 'evoked_data_map_' monk '.mat']);
        resps   = [resps tmp.resps]; %#ok<AGROW>
    end
end


% load the evoked results
load([conf.res_folder 'data' filesep 'evoked_res.mat']);


% shows that the evoked responses tend to reside in the same subspace as
% spanned by the synergies found during natural movement
for i=1:length(resps),
    m1flag(i) = 0;
    if strcmp(resps(i).mapsite(1), 'm1') &  strcmp(resps(i).mapsite(2), 'm1'),
        m1flag(i) = 1; % this is m1
    end
end
projall = zeros(length(resps), 1);

for m = conf.monks

    monk = char(m);

    idx         = strcmp(monk, {resps.monk});
    larger      = find(idx & [resps.field] > 0);
    responses   = vertcat(resps(larger).response);

    proj_res = project(responses, nat_mov_res.(monk).synall_pro, conf.n_boot, conf.noise);
    projall(larger) = proj_res.ratio_dist;
    h = plot_proj(proj_res, m1flag(larger));
    saveas(h, [conf.out 'projection' filesep 'projection_larger0_' monk '.' conf.image_format]);
    close(h);
end

