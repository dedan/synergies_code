
% analysis similarity between results, originally thought to be used only
% to see wether there is a similarity between the evoked responses and the
% synergies found during natural movement. But can also be used, for
% example to look wether pronation or supination data has stronger
% influence on the synergies, found from combined data (pro-sup together)


conf.monks          = {'chalva', 'vega'};
conf.res_folder     = '~/Documents/uni/yifat_lab/results/';
conf.n_boot         = 100;
conf.noise          = 0.5;
conf.image_format   = 'pdf';

resps = struct([]);


% load a result file from the syn_analysis
for monk = conf.monks
    
    if isempty(resps)
        load([conf.res_folder 'data' filesep 'evoked_data_chalva.mat']);     
    else
        tmp     = load([conf.res_folder 'data' filesep 'evoked_data_' char(monk) '.mat']);     
        resps   = [resps tmp.resps]; %#ok<AGROW>
    end
end

% load the nonevoked results
load([conf.res_folder 'data' filesep 'nat_mov_res.mat']);

% load the evoked results
load([conf.res_folder 'data' filesep 'evoked_res.mat']);


%% all evoked on all nonevoked 

% shows that the evoked responses tend to reside in the same subspace as
% spanned by the synergies found during natural movement

for m = conf.monks
    
    monk = char(m);
    
    idx         = strcmp(monk, {resps.monk});
    responses   = vertcat(resps(idx).response);

    proj_res = project(responses, nat_mov_res.(monk).synall, conf.n_boot, conf.noise);
    
    h = plot_proj(proj_res);
    
    saveas(h, [conf.res_folder 'projection' filesep 'projection_' monk '.' conf.image_format]);
    close(h);
    
    disp('principal angles: ');
    
    disp(subspace(evoked_res.(monk).all.nmf', nat_mov_res.(monk).synall'));

end



%%
n = 10000;
r = zeros(1,n);
d = 11;
ds = 3;

for i = 1:n
    r(i) = subspace(rand(d, ds), rand(d,ds));
end
figure
hist(r,100)









    