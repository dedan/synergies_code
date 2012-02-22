

inpath  = '~/projects/yifat_paper/results/data/';
outpath = '~/projects/yifat_paper/results/tuning/';
names   = {'vega', 'darma', 'chalva'};

chan = {'BR-B', 'EDC-E', 'APL-E', 'ECU-E', 'FCR-F', 'APL-E', 'ED45-E', 'ED23-E', ...
            'ECU-E', 'BR-B', 'PL-F', 'FCR-F', 'X', 'X', 'X', 'FDS-F';
        'ECU-E', 'ED45-E', 'EDC-E', 'APL-E', 'ECR-E', 'ED23-E', 'BIC-P', 'BIC-P', ...
            'FDS-F', 'PL-F', 'FCU-F', 'FCR-F', 'PT-F', 'FDP-F', 'TRIC-P', 'BIC-P';
        'FCU-F', 'FDS-F', 'PL-F', 'FCR-F', 'PT-F', 'FDP-F', 'PL-F', 'BIC-P', ...
            'ECU-E', 'EDC-E', 'ED45-E', 'ECR-E', 'ED23-E', 'APL-E', 'ECR-E', 'TRIC-P'};

selection.chalva = 2:15
selection.vega   = 60:80
selection.darma  = 1:10

ex = 1

tuning_x = sort([90 45 0 -45 -90 -135 180 135]*pi/180);
tuning_x = [tuning_x tuning_x(1)];

cur_monk = 4;

% first some statistics
intersection = intersect({chan{1,:}},{chan{2,:}});
intersection = intersect(intersection,{chan{3,:}});
disp('');
disp([num2str(length(tmp)) ' channels in common for all monkeys']);

all_channels = union({chan{1,:}}, {chan{2,:}});
all_channels = union(all_channels, {chan{3,:}});

disp('all channels ever recorded')
disp(all_channels)

for i = 1:length(names)
    disp(names{i});
    disp(intersect({chan{i,:}}, all_channels))
end


% to be sure that the tuning curves are actually consistend of sessions, we
% first pick the first 3 channels and look at their tuning curves over
% sessions

for cur_monk = 1:length(names)

    f = figure('Visible', 'off');
    load([inpath 'all_data_' names{cur_monk}]);
    l = min(20, length(sessions));

    for s = 1:l
    
        used_channels = find(sessions(s).channels);
        for j = 1:3
            h = subplot(3, l, s + l*(j-1));
            bla = used_channels(j);
            tuning_curve = sessions(s).tuning_y{1,bla};
            tuning_curve = [tuning_curve tuning_curve(1)];
            if length(tuning_curve) < 9
                continue
            end
            polar(tuning_x, tuning_curve);
            set(findall(gca, 'Type', 'text'), {'visible'}, {'off'})
        end
    end
    saveas(f, [outpath names{cur_monk} '_overview.png'])
    close(f)
end


% this previous plot might help us to determine sessions with untypical
% tuning to sort them out later. currently we just take all of them!

for cur_monk = 1:length(names)

    load([inpath 'all_data_' names{cur_monk}]);
    c2take   = all(vertcat(sessions.channels));

    mean_pds = circ_mean(vertcat(sessions.pd));
    for c = find(c2take)
        tmp = [];
        
        for s = 1:length(sessions(selection.(names{cur_monk})))
            if length(sessions(s).tuning_y{1, c}) < 8
                continue
            end
            tmp = [tmp; sessions(s).tuning_y{1, c}];
        end
        
        tuning_curve = circ_mean(tmp) ./ sum(circ_mean(tmp));
        tuning_curve = [tuning_curve tuning_curve(1)];

        f = figure('visible', 'off');
        polar(tuning_x, tuning_curve)
        hold on
        polar(mean_pds(c), max(tuning_curve), '*r')
        if length(strmatch(chan{cur_monk, c}, char(chan{cur_monk,:}))) > 1
            saveas(f, [outpath chan{cur_monk, c} '_' names{cur_monk} '_' int2str(ex) '.png']);
            ex = ex + 1;
        else
            saveas(f, [outpath chan{cur_monk, c} '_' names{cur_monk} '.png']);
        end
        close(f)
    end
end


