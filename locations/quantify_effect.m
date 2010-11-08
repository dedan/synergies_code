function   eff = quantify_effect( resp)
% assign colors to responses of different body sites


if ~isempty(findstr( lower(resp), 'finger'))
    eff = 'r';
elseif ~isempty(findstr( lower(resp), 'thumb'))
    eff = 'r';
elseif ~isempty(findstr( lower(resp), 'wrist'))
      eff = 'm';
elseif ~isempty(findstr( lower(resp), 'elbow'))
    eff = 'g';
elseif ~isempty(findstr( lower(resp), 'shoulder'))
    eff = 'b';
elseif ~isempty(findstr( lower(resp), 'sholder'))
    eff = 'b';
else
    disp(['Effect is --> ' resp]);
    eff = 'y';
end
