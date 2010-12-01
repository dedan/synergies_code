function runscript4emgdat(monk)


% path = 'E:\';
path = '/Volumes/LAB/';
config.path     = [path monk filesep];
config.outpath  = [path monk filesep 'EMGdat'];
config.monk     = monk;
config.e2add    = 'e';
config.verbose  = 0;

config.df      = 100;
config.bhv     = 'TO';
config.bck_bhv = 'SC';
config.pre     = -500;
config.post    = 1000;
config.win     = [0 500];           % time relative to torque onset TO
config.bck_win = [-500 0];          % time relative to spatial cue (SC)

% parameters for trial selection (how they are chosen see behav_criteria.m)
% in the data validation folder
config.t_react = [-200, 500];       % reaction time
config.t_move  = [200, 1500];       % movemente time
config.ang_div = 35;                % angular deviation


vdir = dir([config.path 'data' filesep config.monk(1) '*']);
vdir = sortdirs( vdir);

for i= 1:length(vdir),
    
    curdir = char(vdir(i));
    disp(['file ' num2str(i) ' of ' num2str(length(vdir)) ': ' curdir]);
    
    if ~exist([config.outpath filesep 'EMG' curdir '.mat'], 'file'),
        
        [chdata,emgpsth] = MuscleSyn(curdir, config);   %#ok<NASGU>
        if ~isempty(chdata),
            
            
            % give all data the same size (as if recorded from all channels)
            for k = 1:length(chdata)
                for l = 1:length(chdata(k).amp)
                    tmp = zeros(length(chdata(k).channels), size(chdata(k).amp{l},2));
                    tmp(logical(chdata(k).channels),:) = chdata(k).amp{l};
                    chdata(k).amp{l} = tmp;
                    
                    tmp = zeros(length(chdata(k).channels), size(chdata(k).bck_amp{l},2));
                    tmp(logical(chdata(k).channels),:) = chdata(k).bck_amp{l};
                    chdata(k).bck_amp{l} = tmp;
                end
            end
            
            save([config.outpath filesep 'EMG' curdir], 'chdata', 'emgpsth');
            disp(['save EMG' curdir '.mat to ' config.outpath]);
        end
    else
        disp(['EMGdat' filesep 'EMG' curdir '.mat skipped because already exists']);
    end
end

