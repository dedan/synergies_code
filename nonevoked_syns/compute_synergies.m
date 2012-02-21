
% This is a runscript to do all the computationally expensive stuff like
% residual tests. Once this is run, the files which are created can be used
% for further analysis

% call this function for example with: 
% compute_synergies('/Volumes/LAB/', '~/Documents/uni/yifat_lab/results/data/', {'chalva'})

function compute_synergies(data_path, outpath, monks)

addpath('../lib'); %#ok<MCAP>
addpath('../data_validation/');

% general settings
config = struct;
config.Niter_exploration= 50;
config.n_best           = 20;
config.Niter_res_test   = 20;
config.opt              = statset('MaxIter',50);
config.max_channels     = 16;
config.modi             = {'_pro','_sup'};
config.dim              = 3;
config.n_baseline       = 5;

if config.n_best > config.Niter_exploration
    disp('n_best has to be smaller than Niter_exploration');
end

diary([outpath 'log.txt']);

for monk = monks
    config.monk             = char(monk);
    config.emgdat_path      = [data_path config.monk filesep 'EMGdat' filesep];

    tmp_sess = nonevoked_sessions(config); 
    
    disp('computing synergies..');
    sessions = comp_syns(tmp_sess, config);  
    
    disp('estimate matching baselines..');    
    stats    = estimate_baseline(sessions, config); %#ok<NASGU>
    save([outpath 'all_data_' char(monk)], 'sessions', 'stats');
end


diary('off');
disp('finished');



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% estimate the baseline of the matching scores
function stats = estimate_baseline(sessions, conf)

m_scores = NaN(length(sessions), conf.n_baseline);
h_scores = NaN(length(sessions), conf.n_baseline);

c2take   = all(vertcat(sessions.channels));

for i = 1:length(sessions)
    
    data    = sessions(i).mats(1).data_raw(:,c2take);
    hands   = 1;
    
    % if both handpositions available
    if sessions(i).hands > 1 && size(sessions(i).mats(2).data_raw, 1) > 10
        pro     = sessions(i).mats(1).data_raw(:,c2take);
        sup     = sessions(i).mats(2).data_raw(:,c2take);
        data    = vertcat(data, sup);         %#ok<AGROW>
        hands   = 2;
    end
    
    % bootstrap loop
    for j = 1:conf.n_baseline
        r_data          = shuffle_inc(data);
        nmf             = nmf_explore(r_data, conf);
        pca             = pcaica(r_data, conf.dim)';
        [~, ~, scores]  = match_syns(nmf.syns, pca, 1);
        m_scores(i,j)   = mean(scores);
        
        if hands == 2
            r_pro           = shuffle_inc(pro);
            r_sup           = shuffle_inc(sup);
            nmf_pro         = nmf_explore(r_pro, conf);
            nmf_sup         = nmf_explore(r_sup, conf);
            [~, ~, scores]  = match_syns(nmf_pro.syns, nmf_sup.syns, 1);
            h_scores(i,j)   = mean(scores);
        end
    end
end
h_scores = h_scores(~any(isnan(h_scores),2),:);
stats.m_base = mean(m_scores(:));
stats.h_base = mean(h_scores(:));




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function sessions = comp_syns(sessions, conf)

% channels which are available for all sessions
c2take = all(vertcat(sessions.channels));

for i = 1:length(sessions)
    
    disp(['session ' num2str(i)]);

    data    = sessions(i).mats(1).data_raw(:,c2take);
    nmf_res = nmf_explore(data, conf);
    sessions(i).nmf_pro     = nmf_res.syns;
    sessions(i).nmf_pro_std = nmf_res.std;
    sessions(i).pca_pro     = pcaica(data, conf.dim)';
    
    if sessions(i).hands > 1 && size(sessions(i).mats(2).data_raw, 1) > 10
        data    = sessions(i).mats(2).data_raw(:,c2take);
        nmf_res = nmf_explore(data, conf);
        sessions(i).nmf_sup     = nmf_res.syns;
        sessions(i).nmf_sup_std	= nmf_res.std;
        sessions(i).pca_sup     = pcaica(data, conf.dim)';
    else
        sessions(i).hands       = 1;
        sessions(i).nmf_sup     = [];
        sessions(i).nmf_sup_std	= [];
        sessions(i).pca_sup     = [];        
    end
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function res = nonevoked_sessions(config)

% find all sessions and do different residual tests.
% these are later on used to find out only the interesting sessions (not
% already explainable by rank1)

res  = struct;

fdat  = dir([config.emgdat_path 'EMG' config.monk(1) '*.mat']);
rc    = 1;      % result counter

% iterate over all the EMG files
for i= 1:length(fdat)
    
    data = load([config.emgdat_path char(fdat(i).name)]);
    display(['processing file: ' num2str(i) ' of ' num2str(length(fdat))]);
    mats = struct;
    
    % for pronation and supination (middle position skipped because only
    % available for a few sessions
    positions = length(data.chdata);
    if positions > 2
        positions = 2;
    end
    for j=1:positions
        mats(j).data_raw = [];
        
        % for all trials
        for k=1:length(data.chdata(j).amp)
            tmp_amp           = data.chdata(j).amp{k};
            tmp_bck           = data.chdata(j).bck_amp{k};
            
            % normalized by division with acitivity before TO
            mats(j).data_raw  = [mats(j).data_raw; tmp_amp' ./ tmp_bck'];
        end
        
        % remove NaNs from the channels which are not used anyway
        e2take = logical(data.chdata(1).channels);
        mats(j).data_raw(:,~e2take) = 0;
        
        
        % warn if there are still NaNs in the data
        if any(isnan(mats(j).data_raw(:)))
            disp(['nan problem in: ' fdat(i).name]);
            mats(j).data_raw(isnan(mats(j).data_raw)) = 0;
        end
    end
    
    if isempty(mats(1).data_raw)
        disp(['problem with: ' fdat(i).name]);
    else
        
        % put in struct what we got so far
        res(rc).mats     = mats;
        res(rc).name     = char(fdat(i).name);
        res(rc).monk     = config.monk;
        res(rc).hands    = length(data.chdata);
        res(rc).id       = data.chdata.id;
        res(rc).channels = data.chdata.channels;
        res(rc).pd       = vertcat(data.chdata.pd);
        res(rc).tuning_x = vertcat(data.chdata.tuning_x);
        res(rc).tuning_y = vertcat(data.chdata.tuning_y);
        res(rc).p1       = vertcat(data.chdata.p1);
        res(rc).p2       = vertcat(data.chdata.p2);
        res(rc).trials   = data.chdata.trials;
        
        % initialize test results (to make them all same length)
        vars = {'r_nmf', 'std_nmf', 'r_pca', 'r_nmf_s', ...
            'r_nmf_raw', 'std_nmf_raw', 'r_pca_raw', 'r_nmf_s_raw'};
        for k = 1:length(vars)
            for l = 1:positions
                res(rc).([vars{k} config.modi{l}]) = zeros(1,config.max_channels);
            end
        end
        
        
        % explained variance tests
        for j = 1:positions
            c2take  = logical(res(rc).channels);
            [r s]   = test_resid_nmf( mats(j).data_raw(:,c2take), config);
            res(rc).(['r_nmf_raw' config.modi{j}])(1:length(r))    = r;
            res(rc).(['std_nmf_raw' config.modi{j}])(1:length(s))  = s;
            
            r       = test_resid_pcaica( mats(j).data_raw(:,c2take), config);
            res(rc).(['r_pca_raw' config.modi{j}])(1:length(r))    = r;
            
            %repeat for shuffled data as only one random shuffling might not be representative
            tmp = zeros(config.Niter_res_test, min(size(mats(j).data_raw(:,c2take))));
            for k = 1:config.Niter_res_test
                tmp(k,:)    = test_resid_nmf(shuffle_inc(mats(j).data_raw(:,c2take)), config);
            end
            m = mean(tmp);
            res(rc).(['r_nmf_s_raw' config.modi{j}])(1:length(m))  = m;
            
        end
        if positions == 1
            res(rc).('r_nmf_sup')        = 0;
            res(rc).('std_nmf_sup')      = 0;
            res(rc).('r_pca_sup')        = 0;
            res(rc).('r_nmf_s_sup')      = 0;
            res(rc).('r_nmf_raw_sup')    = 0;
            res(rc).('std_nmf_raw_sup')  = 0;
            res(rc).('r_pca_raw_sup')    = 0;
            res(rc).('r_nmf_s_raw_sup')  = 0;
        end
        rc = rc +1;
    end
end




