
% a file for the visual inspection of the emgdata, averaged over trials in
% two different ways and also raw. In the first two plots the y axes
% corresponds two the targets, in the third plot to the trials

% only channels which are available for all sessions of a monkey are used

names = {'chalva', 'vega', 'darma'};
path = '~/projects/yifat_paper/results/';
n = 5

for n = 1:length(names)

    f = figure('Visible', 'off')
    load([path 'data/all_data_' names{n} ]);
    all_chan    = vertcat(sessions.channels);
    c2take      = all(all_chan);

    for i = 1:5
        subplot(2, 5, i);
        imagesc(sessions(i).mats(1).data_raw(:,c2take));
        axis off

        subplot(2, 5, 5  + i);
        dat = sessions(i).mats(1).data_raw(:,c2take);
        imagesc(dat ./ repmat(std(dat), size(dat, 1), 1));
        axis off
        title(num2str(sessions(i).id));

    end
    saveas(f, [path 'validation/data_raw_' names{n} '.png'])
    close(f)
end
