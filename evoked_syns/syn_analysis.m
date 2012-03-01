

%% set confuration

clear

conf = struct;

conf.monks                    = {'chalva', 'vega'};
conf.n_monks                  = length(conf.monks);
conf.modi                     = {'pro', 'sup'};

% general settings
conf.sort_out_zerofield       = 1;        % sort out sessions with 0 muscle response field
conf.skip_window_computation  = 1;
conf.do_resid_test            = 0;
conf.inflag                   = true;     % write info messages
conf.erflag                   = true;     % write error messages
conf.image_format             = 'png';    % format of output images


% emg channels
conf.channels            = [11 12 13 14 21 22 23 24 31 32 33 34 41 42 43 44];

% response
conf.stim_value          = 150;           % look only at stimulations around this value

% matrix factorization
conf.dim                 = 3;             % reduce to this number of synergies
conf.Niter_res_test      = 2;            % number of iterations for the residual test
conf.Niter_exploration   = 3;            % number of iterations for nmf exploration
conf.n_best              = 2;            % take only n_best explorations
conf.opt                 = statset('MaxIter', 50);     % number of early explorations

% clustering
conf.n_cluster = 3;

% normalize over sessions
conf.norm       = true;

% files and folders
conf.result_folder       = '~/projects/yifat_paper/results/';
conf.inpath              = '~/projects/yifat_paper/results/';
conf.cur_res_fold        = [conf.result_folder 'evoked_syns' filesep];
mkdir(conf.cur_res_fold)


% sort sessions according to the handsorted flags
% (see also ../data_validation/sort_sessions)
conf.handsorted          = '';




%% get responses

resps = struct([]);
for monk = conf.monks

    disp(['loading data for: ' char(monk)]);
    if conf.skip_window_computation
        disp('skip window calculation, take data from previous session..');
    else
        disp('calculate responses .. ');
        runscript4evoked('/Volumes/LAB/', monk)
    end

    if isempty(resps)
        load([conf.inpath 'data' filesep 'evoked_data_' char(monk)]);
    else
        tmp     = load([conf.inpath 'data' filesep 'evoked_data_' char(monk)]);
        resps   = [resps tmp.resps]; %#ok<AGROW>
    end
end
clear tmp monk;


%% separate and filter responses

disp(' ');
disp(['total subsessions with stimulation: ' int2str(length(resps))]);

h1 = figure('Visible','off');
h2 = figure('Visible','off');
for m = 1:length(conf.monks)
    monk = char(conf.monks(m));

    % get the indeces of subsessions
    idx.(monk) = strcmp(monk, {resps.monk});
    disp(' ');
    disp([monk ' -- number of sessions: ' num2str(length(find(idx.(monk))))]);


    % sort out the sessions which contain artefacts
    load(['data' filesep conf.handsorted 'sort_' monk]);

    % wow, what a line. I do this weird indexing to make it the same length
    % as flags
    idx.(monk)(idx.(monk)) = idx.(monk)(idx.(monk)) & (flags ~= 2);
    disp(['sorted out ' int2str(length(find(flags == 2))) ...
        ' subsessions because of artefacts (handsorted)']);
    disp([int2str(length(find(idx.(monk)))) ' subsessions remaining']);


    % plot the field size histogram
    figure(h1);
    subplot(length(conf.monks), 1, m)
    [n, xout] = hist([resps(idx.(monk)).field], 0:length(conf.channels));
    bar(xout, n ./ sum(n) * 100);
    title(monk)

    % sort out the sessions with fieldsize 0
    disp(['sorted out ' int2str(length(find([resps(idx.(monk)).field] == 0))) ...
        ' subsessions because of fieldsize was 0']);

    disp(['0-field to non-0-field ratio: ' ...
        num2str( length(find([resps(idx.(monk)).field] == 0)) / length(find(idx.(monk))))]);

    % same indexing magic as above
    idx.(monk)(idx.(monk)) = idx.(monk)(idx.(monk)) & [resps(idx.(monk)).field] ~= 0;

    disp([int2str(length(find(idx.(monk)))) ' subsessions remaining']);

    % check the stimulation amps for the remaining files
    figure(h2);
    amps = abs([resps(idx.(monk)).amp]);
    subplot(length(conf.monks), 1, m)
    hist(amps);
    title(monk);

    % collect responses
    res.(monk).flat = vertcat(resps(idx.(monk)).response);
    res.(monk).field = [resps(idx.(monk)).field];

    % filter responses by significance
    load([conf.inpath 'data' filesep 'nat_mov_res_' monk]);
    c2take      = nat_mov_res.c2take;
    id = find(idx.(monk));
    not_sig_filter = zeros(length(resps(id)), length(conf.channels));
    for i = 1:length(id)
        not_sig_filter(i, resps(id(i)).fields) = 1;
    end
    res.(monk).flat(~logical(not_sig_filter(:, c2take))) = 0;
    % remove "no response lines" (significant muscles were in not in c2take)
    res.(monk).flat = res.(monk).flat(sum(not_sig_filter(:, c2take), 2) ~= 0 ,:);
    res.(monk).field = res.(monk).field(sum(not_sig_filter(:, c2take), 2) ~= 0);

    % normalization
    if conf.norm
        res.(monk).flat = normr(res.(monk).flat);
    end

    disp(['number of recordings: ' num2str(length(find(idx.(monk))))]);
end

saveas(h1, [conf.cur_res_fold  'resp_fields.' conf.image_format]);
close(h1);
saveas(h2, [conf.cur_res_fold  'stim_amp_dist.' conf.image_format]);
close(h2);

disp(' ');
disp('calculated responses');
clear flags monk h1 h2 amps m n xout


%% look at clusters in the responses (imagesc and pca scatter)
colors = {'r', 'g', 'b', 'y'};
colors = {colors{1:conf.n_cluster}}
for i = 1:length(conf.monks)

    monk = conf.monks{i};

    h       = figure('Visible','off');
    dat     = res.(monk).flat;

    subplot(3, 2, 1)
    [cluster_idx c]         = kmeans(dat, conf.n_cluster, 'replicates', 100);
    res.(monk).all.center   = c;
    sizes = histc(cluster_idx, 1:conf.n_cluster);
    [ss, sort_center] = sort(sizes, 'descend');
    bar(ss);
    set(gca,'XTickLabel', {colors{sort_center}});
    title('cluster distribution');

    % plot centers
    subplot(3, 2, 2);
    imagesc(res.(monk).all.center(sort_center, :));
    set(gca,'YTickLabel', {colors{sort_center}});
    set(gca,'XTickLabel', {});
    title([monk ' cluster centers']);

    % project responses on principal components
    [~, loadings] = princomp(dat);
    subplot(3, 2, 3)
    for j = 1:conf.n_cluster

        id = (cluster_idx == j);
        plot(loadings(id,1), loadings(id,2), ['.' colors{j}]);
        hold on
        xlabel('1st PC');
        ylabel('2nd PC');
    end

    % cluster imagesc plot
    subplot(3, 2, 4);
    for j = 1:conf.n_cluster
        tmp(j).c = dat(cluster_idx == j, :);
    end
    agg = [];
    for j = sort_center'
        agg = [agg; tmp(j).c; ones(1, size(dat, 2))];
    end
    for j = 1:length(ss)
        yticks(j) = sum(ss(1:j-1)) + ss(j)/2;
    end
    imagesc(agg);
    set(gca,'YTick', yticks);
    set(gca,'YTickLabel', {colors{sort_center}});
    set(gca,'XTickLabel', {});
    title([monk ' clustered']);

    % kmeans resid test
    subplot(3, 2, 5:6);
    resid = NaN(1,10);
    for j = 1:10
        [~, ~, sumd]  = kmeans(dat, j, 'replicates', 100);
        resid(j)      = sum(sumd);
    end
    plot(resid)
    title('mean field sum')

    saveas(h, [conf.cur_res_fold  'raw_clusters_' monk '.' conf.image_format]);
    close(h);

    f = figure('Visible', 'off');
    for j = 1:length(sort_center)
        subplot(conf.n_cluster, 1, j);
        hist(res.(monk).field(cluster_idx == sort_center(j)), 1:length(conf.channels));
    end
    saveas(h, [conf.cur_res_fold  'cluster_fieldhist_' monk '.' conf.image_format]);
    close(h);

end



%% fieldsize plot

h = figure('Visible','off');
for i = 1:length(conf.monks)

    monk = char(conf.monks(i));
    subplot(length(conf.monks),1,i);

    fields = cell(1, length(conf.channels));
    for j = find(idx.(monk))
        fields{resps(j).field} = [fields{resps(j).field} resps(j).response(resps(j).fields)];
    end

    for j = 1:length(conf.channels);
        hold on
        plot(ones(1, length(fields{j})) * j, fields{j}, '.b');
        plot(j, mean(fields{j}), '.r', 'MarkerSize', 20);
    end
    title(monk);
end
saveas(h, [conf.cur_res_fold  'resp_field_size.' conf.image_format]);
close(h);
clear monk fields h





%% distribution of final data

% plot distribution
h   = figure('Visible','off');
mod = {'pro', 'sup'};

for m = 1:length(conf.monks)
    monk = char(conf.monks{m});


    for k = 1:length(mod)
        subplot(2, length(conf.monks), (m-1) * 2 + k)

        sig_resps   = [];
        unsig_resps = [];
        for i = find(idx.(monk))
            sig                     = false(1, length(resps(i).response));
            sig(resps(i).fields)    = true;
            sig_resps               = [sig_resps resps(i).response(sig)]; %#ok<AGROW>
            unsig_resps             = [unsig_resps resps(i).response(~sig)]; %#ok<AGROW>
        end
        [n_sig, sigx]   = hist(sig_resps, 50);
        [n_unsig, unx]  = hist(unsig_resps, 50);
        n               = sum(n_sig) + sum(n_unsig);
        bar( unx, n_unsig ./ n * 100, 'b')
        hold on
        bar(sigx, n_sig ./ n * 100, 'r')
        hold off
        title([monk ' ' mod{k}]);
    end
end
saveas(h, [conf.cur_res_fold 'response_dist' '.' conf.image_format]);
close(h);
clear h mod monk m k sig_resps unsig_resps sig n_sig sigx n_unsig unx n


%% plot clusters in PD space
for i = 1:conf.n_monks

    % load nonevoked data (PDs are in there)
    load([conf.inpath 'data' filesep 'nat_mov_res_' conf.monks{i}]);

    % roseplot
    h = figure('Visible', 'off');
    center   = res.(conf.monks{i}).all.center;
    pds      = nat_mov_res.pds_all(1,:);    % take only pronation pds

    load('data/channels.mat')
    ch_types = channels.(conf.monks{i}).type(nat_mov_res.c2take);

    plot_rose(h, center, pds, ch_types);
    saveas(h, [conf.cur_res_fold  'center_rose_' conf.monks{i} '.' conf.image_format]);
    close(h);
end
clear i h syn pds rose_agg h_fake j


%% save configuration
evoked_res = res;
save([conf.inpath 'data' filesep 'evoked_res'], 'evoked_res');

disp('finished !');
