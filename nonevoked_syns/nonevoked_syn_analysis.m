
% analysis of the synergist found in a non-stimulation condition
clear

%% add path
addpath('../lib'); %#ok<MCAP>
addpath('../pds'); %#ok<MCAP>


%% set config
run '../config/chalva'


mymap = [linspace(145/255,178/255,32)' linspace(167/255,213/255,32)' linspace(216/255,111/255,32)'];
mymap = [mymap; linspace(178/255,251/255,32)' linspace(213/255,147/255,32)' linspace(111/255,24/255,32)'];
config.mymap = mymap;

mkdir([config.outpath 'sep']);


%% get the data

if config.fetch_data == 1
    all_sessions = nonevoked_sessions(config);
else
    load([config.outpath 'all_nonevoked_sessions']);
end


%% rank1 analysis

% plot the rank1 values for all sessions and the different tests to see
% relation between the different tests

% first store them in a convenient matrix
rank1 = NaN(2,length(all_sessions));
for i = 1:length(all_sessions)
    rank1(1,i)  = all_sessions(i).test.r_nmf_pro(1);
    rank1(2,i)  = all_sessions(i).test.r_pca_pro(1);
end

h = figure('Visible','off');
subplot 211
hold on;

% plot the rank1 values for pca and nmf and also the difference between
% them. furthermore the line at the value at which sessions are chosen as
% significant
plot(rank1','.');
plot(rank1(1,:) - rank1(2,:));
plot(ones(1,length(rank1)) * config.significant);
hold off;
legend('nmf', 'pca', 'dif: nmf - pca');

% also calculate and plot the correlation of the two rank1 values
subplot 212
plot(rank1(1,:),rank1(2,:),'.');
[r p] = corrcoef(rank1');
title(['nmf vs. pcaica, r: ' num2str(r(1,2)) ' - p: ' num2str(p(1,2))]);

saveas(h, [config.outpath  'rank1.' config.image_format]);
close(h);
clear rank1



%% sort out the results which can almost be explained by a rank1 matrix
% NOTE i just realized that I take pro_sup together for this test?
% shouldn't I also separate this?
sig_sessions = struct;
sc           = 1;
for i = 1:length(all_sessions)
    if all_sessions(i).test.r_nmf_pro(1) > config.significant
        sig_sessions(sc).test = all_sessions(i).test;
        sig_sessions(sc).mats = all_sessions(i).mats;
        sig_sessions(sc).info = all_sessions(i).info;
        sc = sc+1;
    end
end
clear sc



%% plot residuals for all sessions in one plot
for i = 1:length(config.names)
    n_max_targets = 8;
    all_resid     = zeros(length(sig_sessions)+2,n_max_targets);
    
    for j = 1:length(sig_sessions)
        tmp_test = sig_sessions(j).test.(['r_nmf' config.names{i}]);
        all_resid(j,1:length(tmp_test)) = tmp_test;
    end
    all_resid(length(sig_sessions)+1,:) = mean(all_resid);
    all_resid(length(sig_sessions)+2,:) = ones(1,size(all_resid,2))*10;
    h = plot_all_resid(all_resid');
    saveas(h, [config.outpath  'all_resid' config.names{i} '.' config.image_format]);
    close(h);

    all_resid     = zeros(length(sig_sessions)+2,n_max_targets);    
    for j = 1:length(sig_sessions)
        tmp_test = sig_sessions(j).test.(['r_nmf_s' config.names{i}]);
        all_resid(j,1:length(tmp_test)) = tmp_test;
    end
    all_resid(length(sig_sessions)+1,:) = mean(all_resid);
    all_resid(length(sig_sessions)+2,:) = ones(1,size(all_resid,2))*10;
    h = plot_all_resid(all_resid');
    saveas(h, [config.outpath  'all_resid_s' config.names{i} '.' config.image_format]);
    close(h);

end

clear all_resid n_max_targets tmp_test


%% plot all at once
if length(sig_sessions(1).mats) > 1
    
    all_resid1 = [];
    all_resid2 = [];
    all_resid3 = [];
    for j = 1:length(sig_sessions)
        all_resid1 = [all_resid1; sig_sessions(j).test.('r_nmf_s_pro')]; %#ok<AGROW>
        all_resid2 = [all_resid2; sig_sessions(j).test.('r_nmf_s_sup')]; %#ok<AGROW>
        all_resid3 = [all_resid3; sig_sessions(j).test.('r_nmf_s')]; %#ok<AGROW>
    end
    
    all_resid1 = [all_resid1; mean(all_resid1)];
    all_resid2 = [all_resid2; mean(all_resid2)];
    all_resid3 = [all_resid3; mean(all_resid3); ones(1,size(all_resid3,2))*10];
    
    h = plot_all_resid_all(all_resid1',all_resid2',all_resid3');
    saveas(h, [config.outpath  'all_resid_all_s' config.names{i} '.' config.image_format]);
    close(h);
end





%% how many synergists do we need to have remaining error smaller then 10 % ?
n_syn_dists = NaN(2, length(sig_sessions));
for i = 1:length(sig_sessions)
    n_syn_dists(1,i)  = find(sig_sessions(i).test.r_nmf_pro < 10, 1);
    n_syn_dists(2,i)  = find(sig_sessions(i).test.r_pca_pro < 10, 1);
end

h = figure('Visible','off');
titles = {'nmf', 'pca'};
for i = 1:size(n_syn_dists,1)
    subplot(2,1,i);
    hist(n_syn_dists(i,:), 3);
    xlabel(titles{i})
end
saveas(h, [config.outpath  'n_syns.' config.image_format]);
close(h);
clear n_syn_dists titles



%% compute the synergies

% from the graphs I chose to take 3 synergies
all_syns = struct;
for i = 1:length(sig_sessions)
    
    tmp_data = sig_sessions(i).mats;
    % compute pcaica synergies and resids
    for j = 1:length(config.names)
        tmp_pca                                             = pcaica(tmp_data(j).data, config.dimensions)';
        sig_sessions(i).res.(['syn_pca' config.names{j}])   = tmp_pca;
        all_syns(i).(['pca' config.names{j}])               = tmp_pca;
        
        
        nmf_res = nmf_explore(tmp_data(j).data, config);

        if config.nmf_stab_plots == 1
            h = figure('Visible','off');
            imagesc(nmf_res.flat);
            title(['standard deviation of group size: ' num2str(nmf_res.std)]);
            saveas(h, [config.outpath  'nmf_expl_stab' int2str(i) '_' int2str(j) '.' config.image_format]);
            close(h);
        end
        
        disp(['standard deviation of group size: ' num2str(nmf_res.std)]);
        
        sig_sessions(i).res.(['syn_nmf' config.names{j}])  = nmf_res.syns;
        all_syns(i).(['nmf' config.names{j}])              = nmf_res.syns;
    end
end
clear tmp_data nmf_res

%% similarity of nmf and pcaica synergies (sep condition)

h = figure('Visible','off');
for i = 1:length(config.names)
    all_nmf = [];
    all_pca = [];
    for j = 1:length(sig_sessions)
        tmp = all_syns(j).(['nmf' config.names{i}])';
        all_nmf = [all_nmf; tmp(:)]; %#ok<AGROW>
        tmp = all_syns(j).(['pca' config.names{i}])';
        all_pca = [all_pca; tmp(:)]; %#ok<AGROW>
    end

    subplot(3,1,i);
    plot(all_nmf, all_pca, '.');
    [r p] = corrcoef([all_nmf all_pca(:)]);
    
    title([config.names{i} ': nmf vs. pca , r: ' num2str(r(1,2)) ' - p: ' num2str(p(1,2))]);
    display([config.names{i} ': nmf vs. pca , r: ' num2str(r(1,2)) ' - p: ' num2str(p(1,2))]);
    xlabel('nmf');
    ylabel('pca');
end
saveas(h, [config.outpath  'nmf_vs_pca_single_sess.' config.image_format]);
close(h);


%% compute synergies, for all sessions at once
all = struct;
fin_syns    = struct;
for i = 1:length(config.names)
    all(i).flat = [];
    for j = 1:length(sig_sessions)
        all(i).flat = [all(i).flat; sig_sessions(j).mats(i).data];
    end
end

for i = 1:length(config.names)
    fin_syns.(['pca_all' config.names{i}]) = pcaica(all(i).flat, config.dimensions)';
    
    nmf_res = nmf_explore(all(i).flat, config);
        
    h = figure('Visible','off');
    imagesc(nmf_res.flat);
    title(['standard deviation of group size: ' num2str(nmf_res.std)]);
    saveas(h, [config.outpath  'nmf_expl_stab_all' int2str(i) '.' config.image_format]);
    close(h);
    
    disp(['standard deviation of group size: ' num2str(nmf_res.std)]);
    
    fin_syns.(['nmf_all' config.names{i}]) = nmf_res.syns;    
end
clear col W0 H0 W H_ score errs bla idx best grouped flat




%% a resid plot for this pooled sessions synergies
for i = 1:length(all)
    resid_res(i,:)  = test_resid_nmf(all(i).flat, config); 
    
    %repeat for shuffled data as only one random shuffling might not be representative
    tmp = zeros(config.Niter_res_test, min(size(all(i).flat)));
    for j = 1:config.Niter_res_test
        tmp(j,:)    = test_resid_nmf(shuffle_inc(all(i).flat), config);
    end
    resid_resr(i,:) = mean(tmp); 
    
end

% add a line to show 10 %
resid_res(i+1,:)    = ones(1,size(resid_res,2)) * 10;
resid_resr(i+1,:)   = ones(1,size(resid_resr,2)) * 10;

% residual test for pcaica
for i = 1:length(all)
    resid_res = [resid_res; test_resid_pcaica(all(i).flat, config)]; %#ok<AGROW>
    
    %repeat for shuffled data as only one random shuffling might not be representative
    tmp = zeros(config.Niter_res_test, min(size(all(i).flat)));
    for j = 1:config.Niter_res_test
        tmp(j,:)    = test_resid_pcaica(shuffle_inc(all(i).flat), config);
    end
    resid_resr      = [resid_resr; mean(tmp)]; %#ok<AGROW>
end

if length(all) > 1
    first = ones(7,1)*100; first(4) = 10;
else
    first = ones(3,1)*100; first(2) = 10;
end
h = createfigure1([first resid_res]', [first resid_resr]',config);
saveas(h, [config.outpath 'resid_all_together.' config.image_format]);
close(h);

clear resid_res resid_resr first tmp



%% consistency over sessions

% nmf seems to be quite stable, nevertheless I take the center of the
% prototypes found in the five runs

for i = 1:length(config.names);
    
    grouped = group(all_syns, ['nmf' config.names{i}]);
    fin_syns.(['nmf' config.names{i}]) = grouped.center;

    
    grouped = group(all_syns, ['pca' config.names{i}]);
    fin_syns.(['pca' config.names{i}]) = grouped.center;
    
    flat = [];
    for j = 1:length(grouped)
        flat = [flat; grouped(j).dat]; %#ok<AGROW>
    end
    
    h = figure('Visible','off');
    subplot(4,1,1:3);
    imagesc(flat);
    colormap(config.mymap);
    axis off
    title(['consistency of synergists over sessions ' config.names{i}]);
    
    subplot(4,1,4);
    imagesc(grouped(1).center);
    colormap(config.mymap);
    axis off
    title(['centers ' config.names{i}]);
    saveas(h, [config.outpath  'syn_consist_sessions_' config.names{i} '.' config.image_format]);
    close(h);
    
    stds(i,:) = grouped(1).idx;   %#ok<SAGROW>
end

h = figure('Visible', 'off');
for i = 1:length(config.names)
    subplot(3,1,i);
    hist(stds(i,:));
    title([config.names{i} ' std of clustering: ' num2str(std(hist(stds(i,:), config.dimensions)))]);
end
saveas(h, [config.outpath  'syn_consist_sessions_std.' config.image_format]);
close(h);


clear flat grouped all_names



%% comparison whether synergies are computed on all sessions together or
% over the single sessions and then grouped and mean taken afterwards

% as it shows a strong similarity of the results, I use the together
% version afterwards because there is less averaging involved
for i = 1:length(config.names)
    h = figure('Visible','off');
    syn_compare(fin_syns.(['nmf' config.names{i}]), fin_syns.(['nmf_all' config.names{i}]), 'single', 'pooled');
    saveas(h, [config.outpath  'all_mean_relation' int2str(i) '.' config.image_format]);
    close(h);
    
end





%% compare the synergists in barplots (pca vs. nmf)

% for single session means
for i = 1:length(config.names)
    
    h = figure('Visible','off');
    syn_compare(fin_syns.(['pca' config.names{i}]), fin_syns.(['nmf' config.names{i}]), 'pca', 'nmf');
    saveas(h, [config.outpath  'best_syn_relation' int2str(i) '.' config.image_format]);
    close(h);
end    

% and for pooled synergies
for i = 1:length(config.names)
    
    h = figure('Visible','off');
    syn_compare(fin_syns.(['pca_all' config.names{i}]), fin_syns.(['nmf_all' config.names{i}]), 'pca', 'nmf');
    saveas(h, [config.outpath  'best_syn_relation_all' int2str(i) '.' config.image_format]);
    close(h);
    
end



%% compare the three best nmf synergists for pronation and supination 
% as nmf and pca showed to lead to very similar results I continued with
% nmf
if length(config.names) > 1
    h = figure('Visible','off');
    syn_compare(fin_syns.nmf_all_pro, fin_syns.nmf_all_sup, 'pronation', 'supination');
    saveas(h, [config.outpath  'between_pos_syn_relation_nmf_all.' config.image_format]);
    close(h);
end

%% look whether I come to the same result if I take pro and sup together
% --> yep, so we take the "together" computation
if length(config.names) > 1
    h = figure('Visible','off');
    syn_compare(fin_syns.nmf_all_pro, fin_syns.nmf_all, 'pronation', 'together');
    saveas(h, [config.outpath  'between_posTog_syn_relation_nmf_all.' config.image_format]);
    close(h);
end
         
%% compare the three best nmf synergists for pronation and supination 
% as nmf and pca showed to lead to very similar results I continued with
% nmf
if length(config.names) > 1
    h = figure('Visible','off');
    syn_compare(fin_syns.nmf_pro, fin_syns.nmf_sup, 'pronation', 'supination');
    saveas(h, [config.outpath  'between_pos_syn_relation_nmf_all_sep.' config.image_format]);
    close(h);
end

%% look whether I come to the same result if I take pro and sup together
% --> yep, so we take the "together" computation
if length(config.names) > 1
    h = figure('Visible','off');
    syn_compare(fin_syns.nmf_pro, fin_syns.nmf, 'pronation', 'together');
    saveas(h, [config.outpath  'between_posTog_syn_relation_nmf_sep.' config.image_format]);
    close(h);
end


%% save the final results

nonevoked_res = struct;
nonevoked_res.all_sessions  = all_sessions;
nonevoked_res.all_syns      = all_syns;
nonevoked_res.config        = config;
nonevoked_res.fin_syns      = fin_syns;
nonevoked_res.sig_sessions  = sig_sessions;
save([config.outpath 'nonevoked_results'], 'nonevoked_res');

%clearvars -except nonevoked_res

