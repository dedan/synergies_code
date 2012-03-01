
% path: path to volume which contains all data
% monks: cell array of names, also cell if only one. example {'chalva'}

function runscript4evoked(path, monks, conf)

if nargin == 2
    disp('no config struct given, standard values used');
    disp('');
    conf.stim_value     = 150;      % look only at stimulations around this value
                                    % set this to [] to use all stimulation values
    conf.window         = [-20 20]; % start and end of average window
    conf.int_window     = [6 20];   % integrate window only in this part (final response)
    conf.inflag         = false;
    conf.errflag        = true;
    conf.channels       = [11 12 13 14 21 22 23 24 31 32 33 34 41 42 43 44];
    disp(conf)
end


conf.inpath         = '~/projects/yifat_paper/results/data/';


for monk = monks
    % load the natural movement results
    load([conf.inpath 'nat_mov_res_' char(monk) '.mat'])

    disp(['calculate responses for ' char(monk) '..']);

    conf.monk               = char(monk);
    conf.dat_folder         = [path char(monk) filesep];
    conf.c2take             = nat_mov_res.c2take;

    % get all subessions in which a stimulation took place
    data = get_all_stimsess(conf, char(monk));

    % filter subessions according to StimAmp
    filtered_data = stimulations_at(data, conf.stim_value);

    resps = responses(filtered_data, conf);

    % add information about session
    for i = 1:length(resps)
        resps(i).c2take = conf.c2take;
        resps(i).monk   = char(monk);
    end
    add = '';
    if isempty(conf.stim_value)
        add = 'all_';
    end
    save([conf.inpath filesep add 'evoked_data_' char(monk)], 'resps');
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function res = get_all_stimsess(config, monk)

% use this function on a folder to get all data about files of sessions in
% which stimulation took place. The resulting list can afterwards be sorted
% by the stimulations_at function

rec = 1;
res = struct('hand', [],'id', {}, 'session', {}, 'subsession', {}, ...
    'file', [], 'amp', [], 'electrode', [], 'location', struct);


dir_list = dir([config.dat_folder 'data' filesep]);
dir_list = sortdirs( dir_list ); % YP sorting the directory list before continuing
for i = 1:length(dir_list)

    filename = dir_list{i};


    % load the file
    load([config.dat_folder 'info_files' filesep filename '_param.mat']);

    for j=1:length(SESSparam.SubSess)

        subs = SESSparam.SubSess(j);

        % which electrodes were used
        % in vega I want to use only stimulation from electrode 2 and 3, 1 is spinal
        used_electrodes = find([DDFparam.Electrode.InUse]);
        if strcmp(monk, 'vega') || strcmp(monk, 'vega_first')
            used_electrodes = used_electrodes(used_electrodes ~= 1);
        end

        for used = used_electrodes
            if isfield(subs.Electrode(used).Stim, 'Flag') && subs.Electrode(used).Stim.Flag
                res(rec).location.depth = subs.Electrode(used).Depth;
                res(rec).location.x     = DDFparam.Electrode(used).X;
                res(rec).location.y     = DDFparam.Electrode(used).Y;
                if isfield(DDFparam.Electrode(used), 'Quad') && isfield(DDFparam, 'Positioner')
                    res(rec).location.quad  = DDFparam.Electrode(used).Quad;
                    res(rec).location.posi  = DDFparam.Positioner;

                end
                res(rec).location.thr   = DDFparam.Electrode(used).Threshold;
                res(rec).location.res   = DDFparam.Electrode(used).Active;
                res(rec).id             = DDFparam.ID;

                res(rec).electrode      = used;
                res(rec).amp            = subs.Electrode(used).Stim.Amp;
                res(rec).hand           = SESSparam.hand;
                res(rec).session        = filename;
                res(rec).subsession     = j;
                res(rec).file           = SESSparam.SubSess(j).Files;

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

% YP: if stim_value is empty all no filtration is applied!

k = 1;
if isempty(stim_value),
    filtered_data = all_stimulations;
    disp('no filteration of stim amp was applied!!');
    return
end

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




