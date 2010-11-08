function filtered = artefact_filter(frq, orig_stimuli, emg_data)

% very simple filter to get rid of the artefact in the stimulation data.
% Was only used in the beginning, now we just skip the part in which the
% artefacdt is located.


% variable to set
window_size = 12;

filter_window_size = window_size * frq;
stimuli = orig_stimuli * frq;

used_for_filter   = floor(length(stimuli)/4);
aggregation       = zeros(used_for_filter, filter_window_size+1);

% learn the filter (just create a average window over some of the stimuli)
for i=1:used_for_filter
   if(stimuli(i) < floor(filter_window_size/2) || floor(stimuli(i)+floor(filter_window_size/2)) > length(emg_data) ), continue; end
   range = floor(stimuli(i)-floor(filter_window_size/2)):floor(stimuli(i)+floor(filter_window_size/2));
   aggregation(i,:) = emg_data(range);
end
filter = mean(aggregation);

% and apply this filter to the data
filtered = emg_data;
for i=1:length(stimuli)
   range = floor(stimuli(i)-floor(filter_window_size/2)):floor(stimuli(i)+floor(filter_window_size/2));
   if(range(1) < 0 || range(end) > length(filtered)), continue; end
   filtered(range) = emg_data(range) - filter;
end

