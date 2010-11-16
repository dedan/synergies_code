function res = nonevoked_sessions(config)

% find all sessions and do different residual tests.
% these are later on used to find out only the interesting sessions (not
% already explainable by rank1)

res  = struct;

fdat  = dir([config.emgdat_path 'EMG' config.monk(1) '*.mat']);
rc    = 1;      % result counter

% iterate over all the EMG files
for i= 1:5 %:length(fdat)
    
    data = load([config.emgdat_path char(fdat(i).name)]);
    display(['processing file: ' num2str(i) ' of ' num2str(length(fdat))]);
    mats = struct;
    
    % for both handpositions
    for j=1:length(data.chdata)
        mats(j).data_raw = [];
        
        for k=1:length(data.chdata.amp)
            tmp_amp           = data.chdata.amp{k};
            tmp_bck           = data.chdata.bck_amp{k};
            
            % normalized by division with acitivity before TO
            mats(j).data_raw  = [mats(j).data_raw; tmp_amp(config.channels2take,:)' ./ tmp_bck(config.channels2take,:)'];
            mats(j).data(k,:) = mean(tmp_amp(config.channels2take,:) ./ tmp_bck(config.channels2take,:),2)';
            mats(j).data_test(k,:) = mean(tmp_amp(config.channels2take,:),2)' ./ mean(tmp_bck(config.channels2take,:),2)';
            

        end
    end
    if length(data.chdata) > 1
        mats(3).data_raw = [mats(1).data_raw; mats(2).data_raw];
        mats(3).data     = [mats(1).data; mats(2).data];
        mats(3).data_test     = [mats(1).data_test; mats(2).data_test];
        
    end
    
    if isempty(mats(1).data)
        disp(['problem with: ' fdat(i).name]);
    else
        
        % put in struct what we got so far
        res(rc).mats           = mats;
        res(rc).info.name      = char(fdat(i).name);
        
        % explained variance tests
        for j = 1:length(mats)
            [r s]                                       = test_resid_nmf( mats(j).data, config);
            res(rc).test.(['r_nmf' config.names{j}])    = r;
            res(rc).test.(['std_nmf' config.names{j}])  = s;
            res(rc).test.(['r_pca' config.names{j}])    = test_resid_pcaica( mats(j).data, config);
            
            %repeat for shuffled data as only one random shuffling might not be representative
            tmp = zeros(config.Niter_res_test, min(size(mats(j).data)));
            for k = 1:config.Niter_res_test
                tmp(k,:)    = test_resid_nmf(shuffle_inc(mats(j).data), config);
            end
            res(rc).test.(['r_nmf_s' config.names{j}]) = mean(tmp);
            
            % and the same for the raw version
            [r s]                                       = test_resid_nmf( mats(j).data_raw, config);
            res(rc).test.(['r_nmf_raw' config.names{j}])    = r;
            res(rc).test.(['std_nmf_raw' config.names{j}])  = s;
            res(rc).test.(['r_pca_raw' config.names{j}])    = test_resid_pcaica( mats(j).data, config);
            
            %repeat for shuffled data as only one random shuffling might not be representative
            tmp = zeros(config.Niter_res_test, min(size(mats(j).data_raw)));
            for k = 1:config.Niter_res_test
                tmp(k,:)    = test_resid_nmf(shuffle_inc(mats(j).data_raw), config);
            end
            res(rc).test.(['r_nmf_s_raw' config.names{j}]) = mean(tmp);

            
        end
        rc = rc +1;
    end
end
all_sessions = res; %#ok<NASGU>

save([config.outpath 'all_nonevoked_sessions'], 'all_sessions');



