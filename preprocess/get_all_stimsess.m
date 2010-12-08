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

