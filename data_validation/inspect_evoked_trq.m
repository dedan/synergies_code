
% look at the stimulation induced torque on a trials by trial basis


clear
close all
addpath('../lib')

%% load the data
path = '/Volumes/LAB/';
monk = 'vega';

load([path 'results' filesep 'data' filesep 'evoked_data_' monk '.mat']);

exp_win     = 50;
replay      = false;
only_sig    = true;

%% explore raw data
for i = 10:length(resps)
    
    for j = 1:length(resps(i).files)
        file = [path monk filesep 'data' filesep resps(i).session filesep ...
            'mat' filesep resps(i).files{j} '_bhv'];
        
        disp(['session: ' num2str(i) ' file: ' num2str(j)]);
        
        stim_check  = whos('-file', file, 'AMstim_on', 'StimTime');
        s_times     = load(file, stim_check.name);
        StimTime    = s_times.(stim_check.name) * 1000;     % convert to ms
        load(file, 'TrqFE', 'TrqRU');
        
        l_win   = (exp_win * 2 +1);
        l       = length(StimTime) * l_win;
        t_fe    = ones(1,l);
        t_ru    = ones(1,l);
        coll_fe = zeros(length(StimTime), l_win);
        coll_ru = zeros(length(StimTime), l_win);
        
        if isempty(StimTime)
            disp('stimtime empty');
            continue;
        end
        
        for k = 1:length(StimTime)
            
            % if not out of bound
            if(floor(StimTime(k) + exp_win) < length(TrqFE) && floor(StimTime(k) - exp_win) > 0);
                
                % compute the range and make it fit
                range = floor(StimTime(k) - exp_win):floor(StimTime(k) + exp_win);
                if size(range,2) > 101
                    range = range(1:101);
                end
                
                % collect all the windows in one long vector
                t_fe((k-1)*l_win +1:k*l_win) = TrqFE(range);
                t_ru((k-1)*l_win +1:k*l_win) = TrqRU(range);
                coll_fe(k,:) = TrqFE(range) - mean(TrqFE(range));
                coll_ru(k,:) = TrqRU(range) - mean(TrqRU(range));
            end
        end
        
        
        % check for significant torque responses by ranksum test
        rank_fe = ranksum(mean(coll_fe(:,30:49),2), mean(coll_fe(:,51:70),2)) < 0.05;
        rank_ru = ranksum(mean(coll_ru(:,30:49),2), mean(coll_ru(:,51:70),2)) < 0.05;
        
        % check for significant torque responses by ttest
        ttest_fe = ttest2(mean(coll_fe(:,30:49),2), mean(coll_fe(:,51:70),2));
        ttest_ru = ttest2(mean(coll_ru(:,30:49),2), mean(coll_ru(:,51:70),2));
        
        if ~(rank_fe || rank_ru || ttest_fe || ttest_ru) && only_sig
            continue;
        end

        figure(1)
        
        if replay
            for k = 1:length(StimTime)
                subplot 311; plot(t_fe(1:1:k*l_win));
                subplot 312; plot(coll_fe(1:k,:)')
                subplot 313; plot(mean(coll_fe(1:k,:)));
                pause(0.0005);
            end
        end
                
        subplot 221
        plot(-exp_win:exp_win, mean(coll_fe));
        title(['mean FE -> ttest: ' num2str(ttest_fe) ' - rank: ' num2str(rank_fe)]);
        set(gca,'ytick',[]) ;
        set(gca,'XTick',-50:20:50)
        
        h1 = subplot(2,2,3);
        plot(t_fe);
        hold on
        plot(exp_win+1:l_win:l, t_fe(exp_win+1:l_win:l), '*r');
        hold off
        set(gca,'ytick',[]) ;
        
        
        subplot 222
        plot(-exp_win:exp_win, mean(coll_ru));
        title(['mean RU -> ttest: ' num2str(ttest_ru) ' - rank: ' num2str(rank_ru)]);
        set(gca,'ytick',[]) ;
        set(gca,'XTick',-50:20:50)
        
        h2 = subplot(2,2,4);
        plot(t_ru);
        hold on
        plot(exp_win+1:l_win:l, t_ru(exp_win+1:l_win:l), '*r');
        hold off
        scrollsubplot(l_win, 1:l, h1, h2);
        set(gca,'ytick',[]) ;
        
        
        pause
    end
end