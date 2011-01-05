
% In this file I started search for something like an "orientation map" on
% the cortex. In the response_map file I already plot a direction which is
% computed by the prefered directions of muscles and a evoked response. Now
% I wanted to see whether the stimulation induced a movement that is
% visible in the torque measurements

clear
close all

%% load the data
path = '/Volumes/LAB/';
monk = 'vega';

load([path 'results' filesep 'data' filesep 'evoked_data_' monk '.mat']);

exp_win     = 50;
sig_win     = 20;
replay      = false;
only_sig    = true;



%% collect the torque measurements around stimulation
res = struct;
for i = 1:length(resps)
    
    res(i).fe = [];
    res(i).ru = [];
    
    for j = 1:length(resps(i).files)
        
        disp(['session: ' num2str(i) ' file: ' num2str(j)]);
        file = [path monk filesep 'data' filesep resps(i).session filesep ...
              'mat' filesep resps(i).files{1} '_bhv'];
        
        
        load(file, 'TrqFE', 'TrqRU');
        stim_check  = whos('-file', file, 'AMstim_on', 'StimTime');
        s_times     = load(file, stim_check.name);
        StimTime    = s_times.(stim_check.name) * 1000;     % convert to ms

        l_win   = (exp_win * 2 +1);
        
        tmp_fe = zeros(length(StimTime), l_win);
        tmp_ru = zeros(length(StimTime), l_win);
        
        for k = 1:length(StimTime)
            
            % if not out of bound
            if(floor(StimTime(k) + exp_win) < length(TrqFE) && floor(StimTime(k) - exp_win) > 0);
                
                % compute the range and make it fit
                range = floor(StimTime(k) - exp_win):floor(StimTime(k) + exp_win);
                if size(range,2) > 101
                    range = range(1:l_win);
                end
                
                % collect windows (mean removed)
                tmp_fe(k,:) = TrqFE(range) - mean(TrqFE(range));
                tmp_ru(k,:) = TrqRU(range) - mean(TrqRU(range));
            end
        end
        
        % connect results from different files
        res(i).fe = [res(i).fe; tmp_fe];
        res(i).ru = [res(i).ru; tmp_ru];
    end
    
    if isempty(res(i).fe)
        disp('alaaaaaaaaaarm');
        continue;
    end
    
    % check for significant torque responses by ranksum and t test
    pre          = mean(res(i).fe(:,exp_win - sig_win:exp_win - 1),2);
    post         = mean(res(i).fe(:,exp_win + 1:exp_win + sig_win),2);
    res(i).sigf  = ranksum(pre, post) < 0.05;
    res(i).sigft = ttest2(pre, post);    
    
    pre          = mean(res(i).ru(:,exp_win - sig_win:exp_win - 1),2);
    post         = mean(res(i).ru(:,exp_win + 1:exp_win + sig_win),2);
    res(i).sigu  = ranksum(pre, post) < 0.05;
    res(i).sigut = ttest2(pre, post);    
end


%% plots

% histogram of significance 
figure(1)
subplot 411
hist([res.sigf])
subplot 412
hist([res.sigu])
subplot 413
hist([res.sigft])
subplot 414
hist([res.sigut])


% by visual inspection I found that the really significant ones were only
% the sessions were t and ranksum test were positive
sig_sessions = find(([res.sigft] & [res.sigf]) | ([res.sigut] & [res.sigu]));
disp(['found ' num2str(length(sig_sessions)) ' significant sessions']);


% plot the significant windows
for i = sig_sessions 

    figure(2)
    % replay the averaging (for debugging purposes)
    if replay
        for k = 1:length(StimTime)
            subplot 311; plot(t_fe(1:1:k*l_win));
            subplot 312; plot(coll_fe(1:k,:)')
            subplot 313; plot(mean(coll_fe(1:k,:)));
            pause(0.0005);
        end
    end
        
    subplot 221
    plot(-exp_win:exp_win, mean(res(i).fe));
    title(['mean FE -> ttest: ' num2str(res(i).sigft) ' - rank: ' num2str(res(i).sigf)]);
    set(gca,'ytick',[]) ;
    set(gca,'XTick',[-exp_win -sig_win 0 sig_win exp_win])
    
    h1      = subplot(2,2,3);
    l       = length(res(i).fe(:));
    flat    = reshape(res(i).fe', [], l);
    plot(flat);
    hold on
    plot(exp_win+1:l_win:l, flat(exp_win+1:l_win:l), '*r');
    hold off
    set(gca,'ytick',[]) ;
    title(['session: ' num2str(i)]);
    
    
    subplot 222
    plot(-exp_win:exp_win, mean(res(i).ru));
    title(['mean FE -> ttest: ' num2str(res(i).sigut) ' - rank: ' num2str(res(i).sigu)]);
    set(gca,'ytick',[]) ;
    set(gca,'XTick',[-exp_win -sig_win 0 sig_win exp_win])
    
    h2      = subplot(2,2,4);
    l       = length(res(i).ru(:));
    flat    = reshape(res(i).ru', [], l);
    plot(flat);
    hold on
    plot(exp_win+1:l_win:l, flat(exp_win+1:l_win:l), '*r');
    hold off
    set(gca,'ytick',[]) ;
    
    scrollsubplot(l_win, 1:l, h1, h2);
    
    
    pause
end









