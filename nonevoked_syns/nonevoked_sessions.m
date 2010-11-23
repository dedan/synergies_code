function res = nonevoked_sessions(config)

% find all sessions and do different residual tests.
% these are later on used to find out only the interesting sessions (not
% already explainable by rank1)

res  = struct;

fdat  = dir([config.emgdat_path 'EMG' config.monk(1) '*.mat']);
rc    = 1;      % result counter

% iterate over all the EMG files
for i= 1:length(fdat)
    
    if exist([config.outpath config.monk(1) num2str(i)], 'file')
        continue;
    end
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
            mats(j).data_raw        = [mats(j).data_raw; tmp_amp' ./ tmp_bck'];
            mats(j).data(k,:)       = mean(tmp_amp ./ tmp_bck,2)';
            mats(j).data_test(k,:)  = mean(tmp_amp,2)' ./ mean(tmp_bck,2)';
        end
        if any(isnan(mats(j).data(:))) || any(isnan(mats(j).data_raw(:)))
            disp(['nan problem in: ' fdat(i).name]);
            mats(j).data(isnan(mats(j).data)) = 0;
            mats(j).data_raw(isnan(mats(j).data_raw)) = 0;
        end
    end
    
    if isempty(mats(1).data)
        disp(['problem with: ' fdat(i).name]);
    else
        
        % put in struct what we got so far
        res(rc).mats     = mats;
        res(rc).name     = char(fdat(i).name);
        res(rc).monk     = config.monk;
        res(rc).id       = data.chdata.id;
        res(rc).hands    = length(data.chdata);
        res(rc).channels = data.chadata(1).channels;
        
        % initialize test results (to make them all same length)
        vars = {'r_nmf', 'std_nmf', 'r_pca', 'r_nmf_s', ...
                'r_nmf_raw', 'std_nmf_raw', 'r_pca_raw', 'r_nmf_s_raw'};
        for k = 1:length(vars)
            for l = 1:positions
                res(rc).([vars{k} config.names{l}]) = zeros(1,config.max_channels);
            end
        end

        
        % explained variance tests 
        for j = 1:positions
            try
            [r s]                                               = test_resid_nmf( mats(j).data, config);
            res(rc).(['r_nmf' config.names{j}])(1:length(r))    = r;
            res(rc).(['std_nmf' config.names{j}])(1:length(s))  = s;
            
            r                                                   = test_resid_pcaica( mats(j).data, config);
            res(rc).(['r_pca' config.names{j}])(1:length(r))    = r;
            
            %repeat for shuffled data as only one random shuffling might not be representative
            tmp = zeros(config.Niter_res_test, min(size(mats(j).data)));
            for k = 1:config.Niter_res_test
                tmp(k,:)    = test_resid_nmf(shuffle_inc(mats(j).data), config);
            end
            m = mean(tmp);
            res(rc).(['r_nmf_s' config.names{j}])(1:length(m))  = m;
            
            % and the same for the raw version
            [r s]                                                   = test_resid_nmf( mats(j).data_raw, config);
            res(rc).(['r_nmf_raw' config.names{j}])(1:length(r))    = r;
            res(rc).(['std_nmf_raw' config.names{j}])(1:length(s))  = s;
            
            r                                                       = test_resid_pcaica( mats(j).data_raw, config);
            res(rc).(['r_pca_raw' config.names{j}])(1:length(r))    = r;
            
            %repeat for shuffled data as only one random shuffling might not be representative
            tmp = zeros(config.Niter_res_test, min(size(mats(j).data_raw)));
            for k = 1:config.Niter_res_test
                tmp(k,:)    = test_resid_nmf(shuffle_inc(mats(j).data_raw), config);
            end
            m = mean(tmp);
            res(rc).(['r_nmf_s_raw' config.names{j}])(1:length(m))  = m;
            
            catch err
                disp(err.message);
                disp('no convergence');
            end 
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
    save([config.outpath config.monk(1) num2str(i)], 'res');
end



