function runscript4emgdat(monk)


path = '/Volumes/LAB/';
config.path     = [path monk filesep];
config.outpath  = [path monk filesep 'EMGdat'];
config.monk     = monk;
config.e2add    = 'e';
config.verbose  = 0;

vdir = dir([config.path 'data' filesep config.monk(1) '*']);
vdir = sortdirs( vdir);

for i= 1 %:length(vdir),

    curdir = char(vdir(i));
    disp(['file ' num2str(i) ' of ' num2str(length(vdir)) ': ' curdir]);
    if ~exist([config.outpath filesep 'EMG' curdir '.mat'], 'file'),
        
        [em2take,e1,e2] = findEMGchannels( curdir, config);
        if ~isempty(em2take),
            [chdata,emgpsth] = MuscleSyn( curdir, e1,e2, em2take, config);  %#ok<NASGU>
            if ~isempty(chdata),
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

