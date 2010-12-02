
% path: path to volume which contains all data
% monks: cell array of names, also cell if only one. example {'chalva'}

function runscript4emgdat(path, monks)

config.e2add    = 'ee';
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

config.also_with_stim   = 0;
config.n_btstrp         = 4000;

for monk = monks
    config.path     = [path char(monk) filesep];
    config.outpath  = [path char(monk) filesep 'EMGdat'];
    config.monk     = char(monk);
    
    vdir = dir([config.path 'data' filesep config.monk(1) '*']);
    vdir = sortdirs( vdir);
    
    for i= 1:length(vdir),
        
        curdir = char(vdir(i));
        disp(['file ' num2str(i) ' of ' num2str(length(vdir)) ': ' curdir]);
        
        if ~exist([config.outpath filesep 'EMG' curdir '.mat'], 'file'),
            
            [chdata, emgpsth] = get_emg4sess(char(vdir(i)), config); %#ok<NASGU>
            if ~isempty(chdata),
                save([config.outpath filesep 'EMG' curdir], 'chdata', 'emgpsth');
                disp(['save EMG' curdir '.mat to ' config.outpath]);
            end
        else
            disp(['EMGdat' filesep 'EMG' curdir '.mat skipped because already exists']);
        end
    end
end

