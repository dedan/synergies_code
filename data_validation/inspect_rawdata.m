
% a file for the visual inspection of the emgdata, averaged over trials in
% two different ways and also raw. In the first two plots the y axes
% corresponds two the targets, in the third plot to the trials

% only channels which are available for all sessions of a monkey are used

path = '~/Documents/uni/yifat_lab/results/natural_mov/';

load([path 'all_data_dummy']);



for i = 1:length(sessions)
    
    disp(['(' num2str(i) ')session id: ' num2str(sessions(i).id) ' - monk: ' sessions(i).monk]);
    
    idx         = strmatch(sessions(i).monk, {sessions.monk});
    all_chan    = vertcat(sessions(idx).channels);
    c2take      = all(all_chan);
    
    subplot 211
    imagesc(sessions(i).mats(1).data(:,c2take));
    title('averaged over trials (per target)');
    
    subplot 212
    imagesc(sessions(i).mats(1).data_raw(:,c2take));
    title('without averaging');

    pause;
end

