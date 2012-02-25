

names = {'chalva', 'vega', 'darma'};
res_path = '~/projects/yifat_paper/results/data_bli/';
inpath = '/Volumes/LAB/';
final = struct

for n = 1:length(names)

    figure
    load([res_path 'all_data_' names{n} ]);
    all_chan    = vertcat(sessions.channels);
    c2take      = all(all_chan);

    res = struct;
    for i = 1:8
        res(i).acu = [];
    end

    fdat  = dir([inpath names{n} filesep 'EMGdat' filesep '*.mat']);
    for i= 1:length(fdat)

        data = load([inpath names{n} filesep 'EMGdat' filesep char(fdat(i).name)]);

        if length(data.chdata(1).amp) < 8
            disp('bla');
            continue
        end

        for j = 1:8
            tmp_amp           = data.chdata(1).amp{j}';
            tmp_bck           = data.chdata(1).bck_amp{j}';
            res(j).acu = [res(j).acu; tmp_amp ./ tmp_bck];
        end
    end

    cum = [];
    for j = 1:8
        cum = [cum; res(j).acu];
    end

    subplot 211
    imagesc(cum)
    title(names{n})

    subplot 212
    imagesc(cum ./ repmat(std(cum), size(cum, 1), 1))

    cum(isnan(cum)) = 0;
    final.(names{n}) = std(cum)


    % saveas(f, [path 'validation/data_raw_' names{n} '.png'])
    % close(f)
end
save(['~/projects/yifat_paper/results/data/norm.mat'], 'final')
