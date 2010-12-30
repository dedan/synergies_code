
% In this file I started search for something like an "orientation map" on
% the cortex. In the response_map file I already plot a direction which is
% computed by the prefered directions of muscles and a evoked response. Now
% I wanted to see whether the stimulation induced a movement that is
% visible in the torque measurements

clear
close all

%% load the data
path = '/Volumes/LAB/';
monk = 'chalva';

load([path 'results' filesep 'data' filesep 'evoked_data_' monk '.mat']);
load([path 'results' filesep 'data' filesep 'nat_mov_res_' monk '.mat']);

exp_win = 50;



%% collect the torque measurements around stimulation
res = struct;
for i = 1:length(resps)
    
    file = [path monk filesep 'data' filesep resps(i).session filesep ...
            'mat' filesep resps(i).files{1} '_bhv'];
        
        
    stim_check  = whos('-file', file, 'AMstim_on', 'StimTime');
    s_times     = load(file, stim_check.name);
    StimTime    = s_times.(stim_check.name);
    load(file, 'TrqFE', 'TrqRU');
    
    if isempty(StimTime)
        res(i).sigf = false;
        res(i).sigu = false;
        res(i).sigft = false;
        res(i).sigut = false;
        continue;
    end
    
    for j = 1:length(StimTime)
        
        % if not out of bound
        if(floor(StimTime(j)*1000 + 100) < length(TrqFE) && floor(StimTime(j)*1000 - 100) > 0);
            
            % compute the range and make it fit
            range = floor(StimTime(j)*1000 - 100):floor(StimTime(j)*1000 + 100);
            if size(range,2) > 201
                range = range(1:201);
            end
            
            figure(1)
            subplot 211
            plot(TrqFE(range));
            subplot 212
            plot(TrqRU(range));
            pause
            
            % collect windows (mean removed)
            res(i).avgf(j,:) = TrqFE(range) - mean(TrqFE(range));
            res(i).avgu(j,:) = TrqRU(range) - mean(TrqRU(range));
            
            % compute direction of movement in window 6 - 20 ms after stim
            len = norm([sum(res(i).avgf(j,106:120)), sum(res(i).avgu(j,106:120))]);
            res(i).x(j)     = sum(res(i).avgf(j,106:120)) / len;
            res(i).y(j)     = sum(res(i).avgu(j,106:120)) / len;
        end
    end
    
    disp(i)
    
    % check for significant torque responses by ranksum test
    p = ranksum(mean(res(i).avgf(:,80:94),2), mean(res(i).avgf(:,106:120),2));
    res(i).sigf = p < 0.05;
    p = ranksum(mean(res(i).avgu(:,80:94),2), mean(res(i).avgu(:,106:120),2));
    res(i).sigu = p < 0.05;

    % check for significant torque responses by ttest
    res(i).sigft = ttest2(mean(res(i).avgf(:,80:94),2), mean(res(i).avgf(:,106:120),2));
    res(i).sigut = ttest2(mean(res(i).avgu(:,80:94),2), mean(res(i).avgu(:,106:120),2));
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


% plot the significant windows
figure(2)

% by visual inspection I found that the really significant ones were only
% the sessions were t and ranksum test were positive
sig_sessions = find(([res.sigft] & [res.sigf]) | ([res.sigut] & [res.sigu]));
disp(['found ' num2str(length(sig_sessions)) ' significant sessions']);

sig_sessions = 1:length(res);
for i = sig_sessions 
    disp(i)
    subplot 311
    plot(-100:100,mean(res(i).avgf));
    subplot 312
    plot(-100:100,mean(res(i).avgu));
    subplot 313
    compass(res(i).x, res(i).y);    
    title(['session: ' num2str(i)]);

    pause;
end



