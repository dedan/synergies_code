function all_stimulations = get_all_stimulations(config)

% use this function on a folder to get all data about files of sessions in
% which stimulation took place. The resulting list can afterwards be sorted
% by the stimulations_at function

rec = 1;
all_stimulations = struct('hand', [],'id', {}, 'session', {}, 'subsession', {}, 'file', [], 'amp', [], 'electrode', [], 'location', struct);


dir_list = dir(config.dat_folder);

for i = 1:length(dir_list)
    
    file = dir_list(i);
    
    % skip hidden and system folders
    if( strcmp(file.name(1), '.')), continue; end;
    
    % load the file
    try
        load([config.dat_folder file.name '/Info/' file.name '_param.mat']);
        
        for j=1:length(SESSparam.SubSess)
            
            subs = SESSparam.SubSess(j);
            ok = 0;
            
            % only if values for stimulation, depth and so on can be found
            if( isfield(subs, 'CT') && ~isempty(subs.CT) && ~isempty(subs.CT.Stim) && ~isempty(subs.CT.Depth) && ~isempty(subs.CT.StimAmp))
                
                all_stimulations(rec).electrode         = 1;
                all_stimulations(rec).amp               = SESSparam.SubSess(j).CT.StimAmp;
                all_stimulations(rec).location.depth    = SESSparam.SubSess(j).CT.Depth;
                if(isfield(DDFparam, 'Cortex'))
                    all_stimulations(rec).location.x      = DDFparam.Cortex.X;
                    all_stimulations(rec).location.y      = DDFparam.Cortex.Y;
                    all_stimulations(rec).location.quad   = DDFparam.Cortex.Quad;
                    ok = 1;
                elseif(isfield(DDFparam, 'Cortex1'))
                    all_stimulations(rec).location.x      = DDFparam.Cortex1.X;
                    all_stimulations(rec).location.y      = DDFparam.Cortex1.Y;
                    all_stimulations(rec).location.quad   = DDFparam.Cortex1.Quad;
                    ok = 1;
                end
            elseif(isfield(subs, 'CT1') && ~isempty(subs.CT1) && ~isempty(subs.CT1.Stim) && ~isempty(subs.CT1.Depth) && ~isempty(subs.CT1.StimAmp))
                
                all_stimulations(rec).electrode       = 2;
                all_stimulations(rec).amp             = SESSparam.SubSess(j).CT1.StimAmp;
                all_stimulations(rec).location.depth  = SESSparam.SubSess(j).CT1.Depth;
                all_stimulations(rec).location.x      = DDFparam.Cortex2.X;
                all_stimulations(rec).location.y      = DDFparam.Cortex2.Y;
                all_stimulations(rec).location.quad   = DDFparam.Cortex2.Quad;
                ok = 1;
            end
            
            % if there is a valid stimulus -> save it
            if(ok)
                all_stimulations(rec).hand            = SESSparam.hand;
                all_stimulations(rec).id              = DDFparam.ID;
                all_stimulations(rec).session         = file.name;
                all_stimulations(rec).subsession      = j;
                all_stimulations(rec).file            = SESSparam.SubSess(j).Files;
                all_stimulations(rec).location.thr    = DDFparam.CTXmap.Thr;
                all_stimulations(rec).location.res    = DDFparam.CTXmap.Resp;
                
                rec = rec +1;
            end
            
        end
        
    catch fehler
        if config.erflag
            switch fehler.identifier
                case {'MATLAB:nonStrucReference'}
                    display([file.name '_param.mat' ' -- SubSess Struct Empty']);
                case {'MATLAB:load:couldNotReadFile'}
                    display(['no Info File for: ' file.name ' available !']);
                otherwise
                    display([file.name ': ' fehler.identifier]);
            end
        end
    end
    
    
end

