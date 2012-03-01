
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

    % load the nonevoked results
    tmp     = load([conf.res_folder 'data' filesep 'nat_mov_res_' monk '.mat']);
    nat_mov_res.(monk) = tmp.nat_mov_res;

    if isempty(resps)
        load([conf.res_folder 'data' filesep 'evoked_data_' monk '.mat']);
        all_r = load([conf.res_folder 'data' filesep 'all_evoked_data_' monk '.mat']);
        all_resps = all_r.resps;
    else
        tmp     = load([conf.res_folder 'data' filesep 'evoked_data_' monk '.mat']);
        resps   = [resps tmp.resps]; %#ok<AGROW>
        all_r   = load([conf.res_folder 'data' filesep 'all_evoked_data_' monk '.mat']);
        all_resps = [all_resps all_r.resps];

    end
end


% load the evoked results
load([conf.res_folder 'data' filesep 'evoked_res.mat']);


%% all evoked on all nonevoked

% shows that the evoked responses tend to reside in the same subspace as
% spanned by the synergies found during natural movement

for m = conf.monks

    monk = char(m);

    idx         = strcmp(monk, {resps.monk});
    larger      = idx & [resps.field] > 0;
    responses   = vertcat(resps(larger).response);
    all_idx         = strcmp(monk, {all_resps.monk});
    all_larger      = all_idx & [all_resps.field] > 0;
    all_responses   = vertcat(all_resps(all_larger).response);

    proj_res = project(responses, all_responses, nat_mov_res.(monk).synall_pro, conf.n_boot, conf.noise);

    h = plot_proj(proj_res);

    saveas(h, [conf.out 'projection' filesep 'projection_larger0_' monk '.' conf.image_format]);
    close(h);


end

