
% path: path to volume which contains all data
% monks: cell array of names, also cell if only one. example {'chalva'}

function runscript4evoked(path, monks)

addpath('../lib');

conf.sampling_rate  = 10000;    % downsampling to this rate (acts as lowpass)
conf.stim_value     = 150;      % look only at stimulations around this value
conf.window         = [-20 20]; % start and end of average window
conf.int_window     = [6 20];   % integrate window only in this part (final response)
conf.inflag         = true;
conf.errflag        = true;
conf.channels       = [11 12 13 14 21 22 23 24 31 32 33 34 41 42 43 44];


conf.inpath         = '~/Documents/uni/yifat_lab/results/data/';

resps = struct([]);

% load the natural movement results
load([conf.inpath 'nat_mov_res.mat'])


for monk = monks
    
    disp(['calculate responses for ' char(monk) '..']);
    
    conf.dat_folder         = [path char(monk) filesep 'data' filesep];
    conf.c2take             = nat_mov_res.(char(monk)).c2take;
    
    % get all subessions in which a stimulation took place
    data = get_all_stimsess(conf);
    
    % filter subessions according to StimAmp
    filtered_data = stimulations_at(data, conf.stim_value);
    
    for i = 1:length(resps)
        resps(i).monk = char(monk);
    end
    
    % collect results
    if isempty(resps)
        resps = responses(filtered_data, conf);
    else
        resps = [resps responses(filtered_data, conf)]; %#ok<AGROW>
    end
end
save([conf.inpath filesep 'evoked_data'], 'resps');




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function res = get_all_stimsess(config)

% use this function on a folder to get all data about files of sessions in
% which stimulation took place. The resulting list can afterwards be sorted
% by the stimulations_at function

rec = 1;
res = struct('hand', [],'id', {}, 'session', {}, 'subsession', {}, ...
    'file', [], 'amp', [], 'electrode', [], 'location', struct);


dir_list = dir(config.dat_folder);

for i = 1:length(dir_list)
    
    file = dir_list(i);
    
    % skip hidden and system folders
    if( strcmp(file.name(1), '.')), continue; end;
    
    % load the file
    load([config.dat_folder file.name filesep 'Info' filesep file.name '_param.mat']);
    
    for j=1:length(SESSparam.SubSess)
        
        subs = SESSparam.SubSess(j);
                
        % which electrodes were used
        % TODO in vega I want to use only stimulation from electrode 2 and 3, 1 is spinal
        used_electrodes = find([DDFparam.Electrode.InUse]);
        for used = used_electrodes
            if subs.Electrode(used).Stim.Flag
                res(rec).electrode      = used;
                res(rec).amp            = subs.Electrode(used).Stim.Amp;
                res(rec).location.depth = subs.Electrode(used).Depth;
                res(rec).location.x     = DDFparam.Electrode(used).X;
                res(rec).location.y     = DDFparam.Electrode(used).Y;
                
                res(rec).hand            = SESSparam.hand;
                res(rec).id              = DDFparam.ID;
                res(rec).session         = file.name;
                res(rec).subsession      = j;
                res(rec).file            = SESSparam.SubSess(j).Files;
                res(rec).location.thr    = DDFparam.Electrode(used).Threshold;
                res(rec).location.res    = DDFparam.Electrode(used).Active;
                
                rec = rec +1;
            end
        end
        
    end
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function filtered_data = stimulations_at(all_stimulations, stim_value)

% this function is used to chose only sessions where the stimulation was
% conducted with a certain value

% counter for result entries and struct declaration
k = 1;
filtered_data = struct('hand', [],'id', {}, 'session', {}, 'subsession', {}, 'file', [], 'amp', [], 'electrode', [], 'location', struct);

% variable initialization for first run of loop
act_depth   = all_stimulations(1).location.depth;
last_depth  = act_depth;
act_file    = all_stimulations(1).session;
last_file   = act_file;
amp_dist    = 10000000;
best_amp_pos = -1;


% loop over all stimulations to filter them around the stim_value
for i = 1:length(all_stimulations)

   session     = all_stimulations(i);
   act_file    = session.session;
   act_depth   = session.location.depth;

   % when the subsession or the depth changes and there
   % is a valid depth -> save it
   if( (~strcmp(last_file,act_file) ||  (last_depth ~= act_depth)) & best_amp_pos ~= -1) %#ok<AND2>
      filtered_data(k) = all_stimulations(best_amp_pos);
      best_amp_pos = i;
      amp_dist    = abs(abs(session.amp) - stim_value);
      k = k +1;

   % nothing changed ? -->  look for a better depth
   else
      if (abs(abs(session.amp) - stim_value) <= amp_dist)
         best_amp_pos = i;
         amp_dist    = abs(abs(session.amp) - stim_value);
      end
   end
   
   %remember the actual file
   last_depth  = act_depth;
   last_file   = act_file;
end

% don't forget the last entry
if(best_amp_pos ~= -1)
   filtered_data(k) = all_stimulations(best_amp_pos);
end




