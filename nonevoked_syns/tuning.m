

inpath  = '~/projects/yifat_paper/results/data/';
outpath = '~/projects/yifat_paper/results/tuning/';
names   = {'vega', 'darma', 'chalva'};
dimensions = 3

selection.chalva = 2:15;
selection.vega   = 60:80;
selection.darma  = 1:10;

ex = 1
tuning_x = sort([90 45 0 -45 -90 -135 180 135]*pi/180);
tuning_x = [tuning_x tuning_x(1)];

% first some statistics
load('data/channels.mat')
disp('');
disp([num2str(length(channels_in_common)) ' channels in common for all monkeys']);
disp('all channels ever recorded')
disp(all_channels)

cic = channels_in_common;
f = figure('Visible', 'off');
for d = 1:dimensions

    for cur_monk = 1:length(names)

        subplot(3, dimensions, (cur_monk-1)*dimensions + d)
        load([inpath 'nat_mov_res_' names{cur_monk}]);

        chan = {channels.(names{cur_monk}).name{nat_mov_res.c2take}};
        cum = zeros(3, length(cic));

        for chan_num = 1:length(cic)

            idx = strmatch(cic{chan_num}, chan);
            for i = 1:length(idx)
                cum(i, chan_num) = nat_mov_res.synnmf_pro(d, idx(i));
            end
        end
        bar(cum')
        title(['synergy #: ' int2str(d) ' - ' names{cur_monk}])
        set(gca,'XTickLabel', char(cic));
    end
end
saveas(f, [outpath 'syn_comparison.png'])


%% synergies in PD space

for cur_monk = 1:length(names)

    load([inpath 'all_data_' names{cur_monk}]);
    load([inpath 'nat_mov_res_' names{cur_monk}]);

    % roseplot
    f = figure('Visible', 'off');
    syn      = nat_mov_res.synnmf_pro;
    pds      = nat_mov_res.pds_all(1,:);

    ch_types = channels.(names{cur_monk}).type(nat_mov_res.c2take);
    plot_rose(f, syn, pds, ch_types);
    saveas(f, [outpath  'syn_rose_pro_' names{cur_monk} '.png']);
    close(f);

    if max([sessions.hands] >1)
        f = figure('Visible', 'off');
        syn      = nat_mov_res.synnmf_sup;
        pds      = nat_mov_res.pds_all(2,:);

        plot_rose(f, syn, pds, ch_types);
        saveas(f, [outpath  'syn_rose_sup_' names{cur_monk} '.png']);
        close(f);
    end
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
    chan_names = channels.(names{cur_monk}).name;

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
        if length(strmatch(chan_names(c), char(chan_names))) > 1
            saveas(f, [outpath chan_names{c} '_' names{cur_monk} '_' int2str(ex) '.png']);
            ex = ex + 1;
        else
            saveas(f, [outpath chan_names{c} '_' names{cur_monk} '.png']);
        end
        close(f)
    end
end


