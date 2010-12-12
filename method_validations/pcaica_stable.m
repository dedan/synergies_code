
% This is a runscript to do all the computationally expensive stuff like
% residual tests. Once this is run, the files which are created can be used
% for further analysis

% call this function for example with: 
% compute_synergies('/Volumes/LAB/', '~/Documents/uni/yifat_lab/results/data/', {'chalva'})

function pcaica_stable()

data_path   = '~/';

addpath('../lib'); %#ok<MCAP>

% general settings
config = struct;
config.monk = 'chalva';
config.max_channels     = 16;
config.modi             = {'_pro','_sup'};
config.dimensions       = 3;


% chalva sessions 54 the connector was switched. 
% channels 1-8 were switched with 9-16
fix_chalva_54(data_path)

config.emgdat_path      = [data_path config.monk filesep 'EMGdat' filesep];

tmp_sess = nonevoked_sessions(config);
sessions = comp_syns(tmp_sess, config);
group_pca = group(sessions, 'pca_pro');
figure
imagesc(vertcat(group_pca.dat));




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function sessions = comp_syns(sessions, conf)

% channels which are available for all sessions
c2take = all(vertcat(sessions.channels));

for i = 1:length(sessions)
    for j = 1:length(sessions(1).hands)
        
        data    = sessions(i).mats(j).data_raw(:,c2take);
        sessions(i).(['pca' conf.modi{j}])          = pcaica(data, conf.dimensions)';
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
            mats(j).data(k,:) = mean(tmp_amp,2)' ./ mean(tmp_bck,2)';
        end
        
        % remove NaNs from the channels which are not used anyway
        e2take = logical(data.chdata(1).channels);
        mats(j).data(:,~e2take) = 0;
        mats(j).data_raw(:,~e2take) = 0;
        
        
        % warn if there are still NaNs in the data
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
        res(rc).hands    = length(data.chdata);
        res(rc).id       = data.chdata.id;
        res(rc).channels = data.chdata.channels;
        res(rc).pd       = data.chdata.pd;
        res(rc).p1       = data.chdata.p1;
        res(rc).p2       = data.chdata.p2;
                
        
        rc = rc +1;
    end
end




