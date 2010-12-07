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
    try
        load([config.dat_folder file.name filesep 'Info' filesep file.name '_param.mat']);
        
        for j=1:length(SESSparam.SubSess)
            
            subs = SESSparam.SubSess(j);
            ok = 0;
            
            % only if values for stimulation, depth and so on can be found
            if( isfield(subs, 'CT') && ~isempty(subs.CT) && ...
             ~isempty(subs.CT.Stim) && ~isempty(subs.CT.Depth) && ...
             ~isempty(subs.CT.StimAmp))
                
                res(rec).electrode         = 1;
                res(rec).amp               = SESSparam.SubSess(j).CT.StimAmp;
                res(rec).location.depth    = SESSparam.SubSess(j).CT.Depth;
                if(isfield(DDFparam, 'Cortex'))
                    res(rec).location.x      = DDFparam.Cortex.X;
                    res(rec).location.y      = DDFparam.Cortex.Y;
                    res(rec).location.quad   = DDFparam.Cortex.Quad;
                    ok = 1;
                elseif(isfield(DDFparam, 'Cortex1'))
                    res(rec).location.x      = DDFparam.Cortex1.X;
                    res(rec).location.y      = DDFparam.Cortex1.Y;
                    res(rec).location.quad   = DDFparam.Cortex1.Quad;
                    ok = 1;
                end
            elseif(isfield(subs, 'CT1') && ~isempty(subs.CT1) && ...
             ~isempty(subs.CT1.Stim) && ~isempty(subs.CT1.Depth) && ...
             ~isempty(subs.CT1.StimAmp))
                
                res(rec).electrode       = 2;
                res(rec).amp             = SESSparam.SubSess(j).CT1.StimAmp;
                res(rec).location.depth  = SESSparam.SubSess(j).CT1.Depth;
                res(rec).location.x      = DDFparam.Cortex2.X;
                res(rec).location.y      = DDFparam.Cortex2.Y;
                res(rec).location.quad   = DDFparam.Cortex2.Quad;
                ok = 1;
            end
            
            % if there is a valid stimulus -> save it
            if(ok)
                res(rec).hand            = SESSparam.hand;
                res(rec).id              = DDFparam.ID;
                res(rec).session         = file.name;
                res(rec).subsession      = j;
                res(rec).file            = SESSparam.SubSess(j).Files;
                res(rec).location.thr    = DDFparam.CTXmap.Thr;
                res(rec).location.res    = DDFparam.CTXmap.Resp;
                
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

