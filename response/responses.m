function resp = responses(filtered_data, config)

% get all the responses (integrated average windows) for all sessions.
% as we have only very little session in the midposition posture, we ignore
% them at all

resp = struct;

j = 1;
for dati = 1:length(filtered_data)
    
    if 1 %(config.inflag == 1)
        disp(' ');
        display(['working on session ' int2str(dati) ' of ' int2str(length(filtered_data))]);
    end
    
    dat = filtered_data(dati);
    
    
    % collect the data for the average window function
    emg_data = [];
    StimTime_agg = [];
    
    % aggregate divided files
    if length(dat.file(1):dat.file(end)) > 1
        disp('have to connect');
    end
    for file_number = dat.file(1):dat.file(end)
        file = [config.dat_folder dat.session '/MAT/' dat.session lead_zeros(file_number,3) int2str(file_number)];
        
        stim_check = whos('-file',[file '_bhv'],'StimTime');
        if(stim_check.size(1) > 4)
            
            
            load([file '_bhv.mat'], 'StimTime');
            load([file '_emg.mat'], 'EMG*');
            first_good = find(config.channels2take,1);
            f_orig = eval(['EMG' int2str(config.channels(first_good)) '_KHz']);
            
            %NOTE beim aneinanderhängen darauf achten zur zweiten stimmtime
            %einen offset in länge der vorigen session draufzurechnen
            
            tmp = NaN(length(find(config.channels2take)), length(eval(['EMG' int2str(config.channels(first_good))])));
            c = 1;
            for i = find(config.channels2take)
                tmp(c,:) = eval(['EMG' int2str(config.channels(i))]);
                c = c+1;
            end
            
            StimTime = StimTime +(size(emg_data,2) / (f_orig*1000));
            StimTime_agg = [StimTime_agg; StimTime]; %#ok<AGROW>
            emg_data = [emg_data tmp]; %#ok<AGROW>
            
            
        end
    end
    
    
    if length(StimTime_agg) > 1
        
        disp('take this session');
        
        resp(j).f_orig      = f_orig;
        resp(j).info        = dat;
        resp(j).hand        = find(dat.hand(:,dat.file(1)));
        resp(j).file        = file;
        resp(j).connected   = length(dat.file(1):dat.file(end)) > 1;
        
        % NOTE continue to skip the non pronation or supination sessions
        if(resp(j).hand ~= 1 && resp(j).hand ~= 2)
            if config.inflag
                display('hand in midposition, skip this session..');
            end
            continue;
        end
        
        wins = get_average_windows(emg_data, StimTime_agg, f_orig, config);
        
        resp(j).x       = wins.x_axis;
        resp(j).windows = wins.windows;
        
        % determine fieldsize (fieldsize is the number of
        % significantly responding channels)
        resp(j).field   = length(find(wins.p < 0.05 & wins.p > 0));
        
        for i = 1:size(wins.windows,1)
            
            % response computation as described by yuval
            resp(j).response_y(i) = (mean(resp(j).windows(i,wins.post_r)) / mean(resp(j).windows(i,wins.pre_r))) * 100;
                        
            % my response computation (integral of average windows and normalization )
            resp(j).response(i)  = sum(abs((resp(j).windows(i,wins.post_r) - wins.m_pre(i)) ./ wins.s_pre(i) ));
            
        end
        
        j = j+1;
    end
end


