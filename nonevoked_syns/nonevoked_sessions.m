function res = nonevoked_sessions(config)

% find all sessions and do different residual tests.
% these are later on used to find out only the interesting sessions (not
% already explainable by rank1)

% only files after the remapping are used (only after and without 250207)
% NOTE I just deleted the files for before 250207

res  = struct;

fdat  = dir([config.emgdat_path 'EMGv*.mat']);
rc    = 1;

for i=1:length(fdat)
    
    data = load([config.emgdat_path char(fdat(i).name)]);
    
    % only when data for different handpositions available
    if length(data.chdata) > 1
        
        display(['processing file: ' num2str(i) ' of ' num2str(length(fdat))]);
        
        % get the data in mat that contains all and in a struct mathand
        % separated by handposition
        mats = struct;
        
        for k=1:2
            mat2add = data.chdata(k).mat';
            if size(data.chdata(k).mat,1) > length(config.channels2take)
                mat2add = mat2add(:,config.channels2take);
            end
            mats(k).data = mat2add;
        end
        mats(3).data = [mats(1).data; mats(2).data];
        
        if ~isempty(mats(1).data)
            
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
            end
            rc = rc +1;
        end
    end
end
all_sessions = res; %#ok<NASGU>

save([config.outpath 'all_nonevoked_sessions'], 'all_sessions');



