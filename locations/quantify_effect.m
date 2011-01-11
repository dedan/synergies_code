% assign colors to responses of different body sites
function   eff = quantify_effect( resp)

% responses were put in the category of the most distal response they were
% containing
fingers     = {'finger', 'fingers', 'thumb', 'fingers (ed-2-3)', ...
               'thumb/finger', 'fingers extension, apl', 'fingers flexion', ...
               'thumb (apl)', 'fingers*,wrist', 'thumbwrist', 'wrist, fingers', ...
               'fingers, wrist', 'fingers, elboe', 'shoulder, fingers'};
wrist       = {'wrist', 'wrist (ecr)', 'wrist*,elbow', 'wrist*, fingers', ...
               'wrist, elbow'};
elbow       = {'elbow', 'elboe', 'wiskers, elbow'};
shoulder    = {'shoulder', 'sholder'};
face_none   = {'face', 'wiskers', 'mouth', 'nose', 'back', ...
               'broke', 'none', 'no response', 'na', 'no stim'};


if is_in(resp, fingers)
    eff = 'r';
elseif is_in(resp, wrist)
    eff = 'g';
elseif is_in(resp, elbow)
    eff = 'b';
elseif is_in(resp, shoulder)
    eff = 'm';
elseif is_in(resp, face_none)
    eff = 'k';
else
    % empty
    eff = 'k';
end
end


function found = is_in(searchstring, c_array)
res     = strfind(c_array, lower(searchstring));
found   = ~isempty([res{:}]);
end
