
function coord = get_cortical_data(ddf, subs, monk)

% which electrodes were used
% in vega I want to use only stimulation from electrode 2 and 3, 1 is spinal
used_electrodes = find([ddf.Electrode.InUse]);
if strcmp(monk, 'vega') || strcmp(monk, 'vega_first')
    used_electrodes = used_electrodes(used_electrodes ~= 1);
end

for used = used_electrodes
    if isfield(subs.Electrode(used).Stim, 'Flag') && subs.Electrode(used).Stim.Flag
        
        if strcmp(monk, 'vega')
            % !!!! X and Y switched for vega
            coord.x         = ddf.Electrode(used).Y;
            coord.y         = ddf.Electrode(used).X;
            coord.electrode = used;
            coord.posi      = ddf.Positioner;
            coord.qd        = ddf.Electrode(used).Quad;
        else
            coord.x         = ddf.Electrode(used).X;
            coord.y         = ddf.Electrode(used).Y;
        end
    end
end