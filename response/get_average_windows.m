function res = get_average_windows(emg_data, StimTime, f_orig, config)

% the function computes average windows for all emg channels around
% stimulations. Configurations are set in the config struct parameter.
% The EMG data is first rectified and dc-value for each channel
% is removed. Then the windows are computed and afterwards dc-removal and a
% normalization by standard-deviation is done for each channel seperately.
% (if set in config)
% this was the old procedure, now we don't normalize, but take the
% percentage between average pre stimulus activation and average post
% stimulus activation

res = struct;

% variable calculation and initialization
sample_dist_ms       = (1 / config.sampling_rate) * 1000;
res.x_axis           = config.window(1)  : sample_dist_ms : config.window(2);
res.windows          = zeros(size(emg_data,1), length(res.x_axis));


try
    
    
    % the values have to be computed to fit the original sampling rate
    StimTime_ms = StimTime * 1000;
    stim_time_up  = StimTime_ms * f_orig;
    win_start_up  = config.window(1) * f_orig;
    win_end_up    = config.window(2)   * f_orig;
    window_size_up = abs(win_end_up - win_start_up)+1;
    aggregation_matrix   = zeros(length(StimTime), window_size_up);
    
    
    if config.inflag
        display('calculating the averages for the timewindows..');
    end
    
    % computation
    for i = 1:size(emg_data,1);
        
        % make tmp variable for channel
        cur_chan = emg_data(i,:);
        
        % DC removal, rectification
        cur_chan = abs(cur_chan - mean(cur_chan));
        
        % apply the artefact filter if flag is set
        if(config.apply_filter)
            cur_chan = artefact_filter(f_orig, StimTime_ms, cur_chan);
        end
        
        % create window for all stimulations
        for j = 1:length(StimTime)
            
            % if window valid, aggregate the window
            if(stim_time_up(j) + win_start_up > 0 && stim_time_up(j) + win_end_up < length(cur_chan))
                range = floor(stim_time_up(j)+win_start_up) : floor(stim_time_up(j)+win_end_up);
                aggregation_matrix(j,:) = cur_chan( range(1:window_size_up));
            end
        end
        
        % check wether responses is of this channel is signifikant like
        % discribed in yuval's short term report
        res.pre_r   = 1:floor(size(aggregation_matrix,2)/2) - config.int_window(1) * f_orig;
        res.post_r  = floor(size(aggregation_matrix,2)/2) + config.int_window(1) * f_orig:size(aggregation_matrix,2);
        res.p(i)    = ranksum(mean(aggregation_matrix(:,res.pre_r),2), mean(aggregation_matrix(:,res.post_r),2));
        
        % compute the mean in the pre stimulation part of the window to
        % remove the individual dc component of a channel. Also std for
        % this part is computed for normalization
        pre_stim = aggregation_matrix(:,res.pre_r);
        res.m_pre(i) = mean(pre_stim(:));
        res.s_pre(i) = std(pre_stim(:));
        
        % average over responses
        res.windows(i,:) = mean(aggregation_matrix);
        
        
    end
    
    if config.inflag
        disp(' ');
        disp('averages computed !');
    end
    
    % error handling
catch fehler
    if config.erflag
        disp('average windows error handling: ');
        display(['Error: ' fehler.identifier]);
        error('error in average window');
    end
end


