
% look at the stimulation induced torque on a trials by trial basis


clear
close all

%% load the data
path = '/Volumes/LAB/';
monk = 'chalva';

load([path 'results' filesep 'data' filesep 'evoked_data_' monk '.mat']);

exp_win = 50;

%% explore raw data
for i = 1:length(resps)
    
    file = [path monk filesep 'data' filesep resps(i).session filesep ...
            'mat' filesep resps(i).files{1} '_bhv'];
        
        
    stim_check  = whos('-file', file, 'AMstim_on', 'StimTime');
    s_times     = load(file, stim_check.name);
    StimTime    = s_times.(stim_check.name) * 1000;     % convert to ms
    load(file, 'TrqFE', 'TrqRU');
       
    l_win   = (exp_win * 2 +1);
    l       = length(StimTime) * l_win;
    t_fe    = ones(1,l);
    t_ru    = ones(1,l);
    
    for j = 1:length(StimTime)
        
        % if not out of bound
        if(floor(StimTime(j) + exp_win) < length(TrqFE) && floor(StimTime(j) - exp_win) > 0);
            
            % compute the range and make it fit
            range = floor(StimTime(j) - exp_win):floor(StimTime(j) + exp_win);
            if size(range,2) > 101
                range = range(1:101);
            end
            
            % collect all the windows in one long vector
            t_fe((j-1)*l_win +1:j*l_win) = TrqFE(range);
            t_ru((j-1)*l_win +1:j*l_win) = TrqRU(range);
        end
    end
    
    % plot
    figure(1)
    h1 = subplot(2,1,1);
    plot(t_fe);
    hold on
    plot(exp_win+1:l_win:l, t_fe(exp_win+1:l_win:l), '*r');
    hold off
    
    h2 = subplot(2,1,2);
    plot(t_ru);
    hold on
    plot(exp_win+1:l_win:l, t_ru(exp_win+1:l_win:l), '*r');
    hold off   
    scrollsubplot(l_win, 1:l, h1, h2);
    pause

end