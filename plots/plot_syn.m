function plot_syn(f, syn, monkey_channels)

    extensors    = {'ECR-E'  'EDC-E'  'ECU-E'  'ED23-E' 'ED45-E'};
    flexors      = {'FCR-F'  'PL-F'   'FCU-F'  'FDS-F'  'FDP-F'};
    others       = {'APL-E'  'TRIC-P' 'PT-F'   'BIC-P'  'BR-B'};
    all_channels = [extensors others flexors];
    ids = {1:5, 6:10, 11:15};
    colors = 'rgb';

    cum = zeros(3, length(all_channels));

    for chan_num = 1:length(all_channels)

        idx = strmatch(all_channels{chan_num}, monkey_channels);
        for i = 1:length(idx)
            cum(i, chan_num) = syn(idx(i));
        end
    end
    bar(cum')

    % rotate tick labels
    h = gca;
    set(h,'XTickLabel',[]);
    b=get(h,'XTick');
    c=get(h,'YTick');
    for i = 1:length(colors)
        tmp_b = b(ids{i});
        tmp_c = repmat(c(1)-.1*(c(2)-c(1)), length(tmp_b),1);
        text(tmp_b, tmp_c, char(all_channels{ids{i}}), ...
            'HorizontalAlignment', 'right', ...
            'rotation', 45, ...
            'color', colors(i));
    end



