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


