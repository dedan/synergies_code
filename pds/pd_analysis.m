
%%
% with a roseplot, where muscle activation patterns are plotted in the
% space of preferred directions of their corresponding muscles, it can be
% seen wether the found synergies are also biologically reasonable.
% a muscle synergy which leads to a activation of muscles, all pointing in
% different directions does not make sense..


%% load data from previous computations

addpath('../lib'); %#ok<MCAP>


% configuration
config = struct;
config.dimensions       = 3;
config.channels2take    = find(logical([0  1  1  1  1  1  1  1  1  1  1  1  0  0  0  0 ]));
config.n_used_channels  = length(config.channels2take);
config.emgdat_path      = '/Volumes/DEDAN_DISK/vega/EMGdat/';
config.outpath          = '~/Documents/uni/yifat_lab/results/pds/';
config.image_format     = 'pdf';
config.res_folder       = '~/Documents/uni/yifat_lab/results/';
config.pd_folder        = '/Volumes/DEDAN_DISK/vega/pd_files/';
config.names            = {'evo', 'non_mean','non_all'};
config.postures         = {'pro', 'sup'};
config.n_bootstrap      = 100;

% load results of previous computations
load([config.res_folder 'nonevoked_syns' filesep 'nonevoked_results']);
all_syns = nonevoked_res.all_syns;
fin_syns = nonevoked_res.fin_syns;
load([config.res_folder 'data' filesep 'latest_evoked']);

all_synergies = struct;

% synergies from stimulation
all_synergies.evo_pro = fin_res(1).nmf_syns ./ max(max(fin_res(1).nmf_syns));
all_synergies.evo_sup = fin_res(2).nmf_syns ./ max(max(fin_res(2).nmf_syns));

% nonevoked synergies
all_synergies.non_mean_pro = fin_syns.nmf_pro ./ max(max(fin_syns.nmf_pro));
all_synergies.non_mean_sup = fin_syns.nmf_sup ./ max(max(fin_syns.nmf_sup));
all_synergies.non_all_pro = fin_syns.nmf_all_pro ./ max(max(fin_syns.nmf_all_pro));
all_synergies.non_all_sup = fin_syns.nmf_all_sup ./ max(max(fin_syns.nmf_all_sup));


% get pds
pds = get_pds(config);

%%

rose(pds(1,:),360)



%% create roseplots of synergists



for i = 1:length(config.names)
    
    
    % rose plots
    h = figure('Visible','off');
    for j = 1:config.dimensions
        
        % roseplot pronation
        subplot(3,2,((j-1)*2)+1);
        rose_agg = [];
        syn = all_synergies.([config.names{i} '_pro']);
        for k = 1:config.n_used_channels
            rose_agg = [rose_agg ones(1, floor(syn(j,k) * 100)) * pds(1,k)]; %#ok<AGROW>
        end
        h_fake = rose(ones(1,100));
        hold on;
        rose(rose_agg);
        set(h_fake, 'Visible', 'Off');
        
        % roseplot supination
        subplot(3,2,((j-1)*2)+2);
        rose_agg = [];
        syn = all_synergies.([config.names{i} '_pro']);
        for k = 1:config.n_used_channels
            rose_agg = [rose_agg ones(1, floor(syn(j,k) * 100)) * pds(2,k)]; %#ok<AGROW>
        end
        h_fake = rose(ones(1,100));
        hold on;
        rose(rose_agg);
        set(h_fake, 'Visible', 'Off');
    end
    saveas(h, [config.outpath 'finsyn_rose_' config.names{i} '.' config.image_format]);
    close(h);
    
end





%% test

% da alles nicht signifkant und scheisse ist was da unten noch kommt, schaue ich mal ob diese pd und
% cstd berechnung denn wenigstens für die original aktivierungsmuster was ordentliches rausbekommt.


% pronation supination
for i = 1:2%length(config.postures)
    
    c = 0;
    
    % get distribution of muscle activations for comparison (bootstrap)
    flat = [];
    for j = 1:length(nonevoked_res.sig_sessions)
        flat = [flat; nonevoked_res.sig_sessions(j).mats(i).data(:)]; %#ok<AGROW>
    end
    flat = flat(:);
    
    
    % for all sessions
    for j = 1:length(nonevoked_res.sig_sessions)
        
        % for all synergies
        for k = 1:config.dimensions
            
            c = c+1;
            
            % bootstrap
            for boot=1:config.n_bootstrap;
                perms = randperm(size(flat,1));
                rand_act = flat(perms(1:config.n_used_channels));
                rand_act(8) = 0;                                                    % discard channel 9
                boot_dist(c,boot) = circ_std( pds(i, rand_act > median(rand_act))); %#ok<AGROW>
            end
            
            
            syn = nonevoked_res.sig_sessions(j).mats(i).data(floor(rand(1)*7)+1,:);
            syn(8) = 0;                                                             % discard channel 9
            cstd(c) = circ_std( pds(i, syn > median(syn))); %#ok<AGROW>
        end
    end
    handle = figure('Visible','off');
    [a b] = hist(boot_dist(:));
    a = a./config.n_bootstrap;
    bar(b,a,'b');
    hold on;
    [a b] = hist(cstd);
    bar(b,a,'r');
    hold off;
    
    [p h] = ranksum(cstd,boot_dist(:,1));
    title(['significant difference: ' int2str(h) ' p-value: ' num2str(p)]);
    saveas(handle, [config.outpath 'circ_std_orig_dat_test_' config.postures{i} '.' config.image_format]);
    close(handle);
    clear cstd boot_dist;
end



%% computation of circular std deviation, for synergies of single sessions
% channel 8 (of used channels) discarded as no reliable pd information available (see
% feather plot of pd consistency)



% pronation supination
for i = 1:2%length(config.postures)
    
    c = 0;
    
    % get distribution of muscle activations for comparison (bootstrap)
    flat = [];
    for j = 1:length(all_syns)
        flat = [flat; all_syns(j).(['nmf_' config.postures{i}])]; %#ok<AGROW>
    end
    flat = flat(:);
    
    
    % for all sessions
    for j = 1:length(all_syns)
        
        % for all synergies
        for k = 1:config.dimensions
            
            c = c+1;
            
            % bootstrap
            for boot=1:config.n_bootstrap;
                perms = randperm(size(flat,1));
                rand_act = flat(perms(1:config.n_used_channels));
                rand_act(8) = 0;                                                    % discard channel 9
                boot_dist(c,boot) = circ_std( pds(i, rand_act > median(rand_act))); %#ok<AGROW>
            end
            
            
            syn = all_syns(j).(['nmf_' config.postures{i}])(k,:);
            syn(8) = 0;                                                             % discard channel 9
            cstd(c) = circ_std( pds(i, syn > median(syn))); %#ok<AGROW>
        end
    end
    handle = figure('Visible','off');
    [a b] = hist(boot_dist(:));
    a = a./config.n_bootstrap;
    bar(b,a,'b');
    hold on;
    [a b] = hist(cstd);
    bar(b,a,'r');
    hold off;
    
    [p h] = ranksum(cstd,boot_dist(:,1));
    title(['significant difference: ' int2str(h) ' p-value: ' num2str(p)]);
    saveas(handle, [config.outpath 'circ_std_persessions_' config.postures{i} '.' config.image_format]);
    close(handle);
    clear cstd boot_dist;
end



%% for the final synergies


% different synergies
for s = 1:length(config.names);
    
    c = 0;
    % pronation supination
    for i = 1:length(config.postures)
        
        
        % get distribution of muscle activations for comparison (bootstrap)
        flat = [];
        for j = 1:length(all_syns)
            flat = [flat; all_syns(j).(['nmf_' config.postures{i}])]; %#ok<AGROW>
        end
        flat = flat(:);
        
        % for all synergies
        for k = 1:config.dimensions
            
            c = c+1;
            
            % bootstrap
            for boot=1:config.n_bootstrap;
                perms = randperm(size(flat,1));
                rand_act = flat(perms(1:config.n_used_channels));
                rand_act(8) = 0;
                boot_dist(boot) = circ_std( pds(i, rand_act > median(rand_act))); %#ok<AGROW>
            end
            
            
            syn = all_synergies.([config.names{s} '_' config.postures{i}])(k,:);
            syn(8) = 0;                                                             % discard channel 9
            cstd(c) = circ_std( pds(i, syn > median(syn))); %#ok<AGROW>
        end
        
    end
    handle = figure('Visible','off');
    
    [a b] = hist(boot_dist(:));
    a = a./config.n_bootstrap;
    bar(b,a,'b');
    hold on;
    [a b] = hist(cstd);
    bar(b,a,'r');
    hold off;
    
    [p h] = ranksum(cstd,boot_dist(:,1));    
    
    
    title(['significant difference: ' int2str(h) ' p-value: ' num2str(p)]);
    saveas(handle, [config.outpath 'circ_std_syns_' config.names{s} '.' config.image_format]);
    close(handle);
    clear cstd boot_dist;
    
end

