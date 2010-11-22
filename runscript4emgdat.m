function runscript4emgdat(monk)


path = 'E:\';
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
config.win     = [0 500];  % time relative to torque onset TO
config.bck_win = [-500 0]; % time relative to spatial cue (SC)

all = [11 12 13 14 21 22 23 24 31 32 33 34 41 42 43 44];

vdir = dir([config.path 'data' filesep config.monk(1) '*']);
vdir = sortdirs( vdir);

for i= 1:length(vdir),

    curdir = char(vdir(i));
    disp(['file ' num2str(i) ' of ' num2str(length(vdir)) ': ' curdir]);
    if ~exist([config.outpath filesep 'EMG' curdir '.mat'], 'file'),
        
        [em2take,e1,e2] = findEMGchannels( curdir, config);
        
        channels = zeros(1,length(all));
        for j = 1:length(em2take)
            channels(all == em2take(j)) = 1;
        end           
        
        if ~isempty(em2take),
            [chdata,emgpsth] = MuscleSyn( curdir, e1,e2, em2take, config);  %#ok<NASGU>
            if ~isempty(chdata),
                chdata.channels = channels;
                save([config.outpath filesep 'EMG' curdir], 'chdata', 'emgpsth');
                disp(['save EMG' curdir '.mat to ' config.outpath]);
            else
                disp('No available data');
            end
        end
    else
        disp(['EMGdat' filesep 'EMG' curdir '.mat skipped because already exists']);
    end
end

