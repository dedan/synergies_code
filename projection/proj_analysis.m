
% analysis similarity between results, originally thought to be used only
% to see wether there is a similarity between the evoked responses and the
% synergies found during natural movement. But can also be used, for
% example to look wether pronation or supination data has stronger
% influence on the synergies, found from combined data (pro-sup together)


res_folder  = '~/Documents/uni/yifat_lab/results/';
n_boot      = 10000;
noise       = 0.5;

% load a result file from the syn_analysis
load('~/Documents/uni/yifat_lab/results/data/040909_1529.mat');     
% 270109_1829

% load the nonevoked results
load([res_folder 'nonevoked_syns/nonevoked_results.mat']);
close all;



%% all evoked on all nonevoked 

% shows that the evoked responses tend to reside in the same subspace as
% spanned by the synergies found during natural movement


proj_res = project(fin_res(3).dat', nonevoked_res.fin_syns.nmf_all', n_boot, noise);
figure
plot_proj(proj_res);

    