

conf.inpath   = '~/projects/yifat_paper/results/data/';
outpath = '~/projects/yifat_paper/results/tuning/';

names = {'vega_first', 'vega_later', 'darma', 'chalva'};

chan = {'FCU-F', 'X', 'PL-F', 'FCR-F', 'BR-B', 'X', 'FDS-F', 'FDP-F', ...
            'X', 'EDC-E', 'X', 'APL-E', 'EDC-E', 'ECR-E', 'ED45-E', 'ECU-E';
        'BR-B', 'EDC-E', 'APL-E', 'ECU-E', 'FCR-F', 'APL-E', 'ED45-E', 'ED23-E', ...
            'ECU-E', 'BR-B', 'PL-F', 'FCR-F', 'X', 'X', 'X', 'FDS-F';
        'ECU-E', 'ED45-E', 'EDC-E', 'APL-E', 'ECR-E', 'ED23-E', 'BIC-P', 'BIC-P', ...
            'FDS-F', 'PL-F', 'FCU-F', 'FCR-F', 'PT-F', 'FDP-F', 'TRIC-P', 'BIC-P';
        'FCU-F', 'FDS-F', 'PL-F', 'FCR-F', 'PT-F', 'FDP-F', 'PL-F', 'BIC-P', ...
            'ECU-E', 'EDC-E', 'ED45-E', 'ECR-E', 'ED23-E', 'APL-E', 'ECR-E', 'TRIC-P'};

cur_monk = 4;

% to be sure that the tuning curves are actually consistend of sessions, we
% first pick the first 3 channels and look at their tuning curves over
% sessions
        
figure(1)
load([conf.inpath 'all_data_' names{cur_monk}]);

for s = 1:length(sessions)
    sess = sessions(s);

    used_channels = find(sess.channels);
    for j = 1:3
        h = subplot(3, length(sessions), s + length(sessions)*(j-1));
        bla = used_channels(j);
        polar(sess.tuning_x{1,bla}, sess.tuning_y{1,bla});
        set(findall(gca, 'Type', 'text'), {'visible'}, {'off'})
        p = get(h, 'pos');
        p(4) = p(4) + 0.1;
        p(3) = p(3) + 0.01;
        set(h, 'pos', p);
    end
end

% from this plot the session ids we want to average over, e.g. for chalva
% we decided to average over all but the first session

sessions = sessions(2:end);

for c = find(sessions(1).channels)
    tmp = [];
    
    for s = 1:length(sessions)
        tmp = [tmp; sessions(s).tuning_y{1, c}];
    end
    tuning_curve = mean(tmp);
    f = figure('visible', 'off');
    polar(sessions(1).tuning_x{1,1}, tuning_curve ./ sum(tuning_curve))
    title(chan{cur_monk, c})
    saveas(f, [outpath names{cur_monk} '_' int2str(c) '.png'])
    close(f)
end
        
        
    
