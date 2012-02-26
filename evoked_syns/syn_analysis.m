

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
conf.image_format             = 'pdf';    % format of output images


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
    load([conf.inpath 'data' filesep conf.handsorted 'sort_' monk]);

    % wow, what a line. I do this weird indexing to make it the same length
    % as flags
    idx.(monk)(idx.(monk)) = idx.(monk)(idx.(monk)) & (flags ~= 2);
    disp(['sorted out ' int2str(length(find(flags == 2))) ...
        ' subsessions because of artefacts (handsorted)']);
    disp([int2str(length(find(idx.(monk)))) ' subsessions remaining']);


    % plot the field size histogram
    figure(h1)
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


    % sort the resulting vectors according to the hand position of the subsession
    idx.pro = ([resps.hand] == 1);
    idx.sup = ([resps.hand] == 2);

    if conf.norm
        res.(monk).pro.flat = normr(vertcat(resps(idx.(monk) & idx.pro).response));
        if ~isempty(vertcat(resps(idx.(monk) & idx.sup).response))
            res.(monk).sup.flat = normr(vertcat(resps(idx.(monk) & idx.sup).response));
        else
            res.(monk).sup.flat = [];
        end
    else
        res.(monk).pro.flat = vertcat(resps(idx.(monk) & idx.pro).response);
        if ~isempty(vertcat(resps(idx.(monk) & idx.sup).response))
            res.(monk).sup.flat = vertcat(resps(idx.(monk) & idx.sup).response);
        else
            res.(monk).sup.flat = [];
        end
    end
    res.(monk).all.flat = vertcat(res.(monk).pro.flat, res.(monk).sup.flat);

    disp(['pronation - number of recordings: ' num2str(length(find(idx.(monk) & idx.pro)))]);
    disp(['supination - number of recordings: ' num2str(length(find(idx.(monk) & idx.sup)))]);

end

saveas(h1, [conf.cur_res_fold  'resp_fields.' conf.image_format]);
close(h1);
saveas(h2, [conf.cur_res_fold  'stim_amp_dist.' conf.image_format]);
close(h2);

disp(' ');
disp('calculated and separated responses');
clear flags monk h1 h2 amps m n xout





%% look at clusters in the responses (imagesc and pca scatter)
colors           = {'r', 'g', 'b', 'y'};

% without loop stupid code separation of monkeys because so different

%% for chalva
monk    = 'chalva';
h       = figure('Visible','off');
dat     = res.(monk).pro.flat;

subplot(3, 3, 1)
[cluster_idx c]         = kmeans(dat, conf.n_cluster, 'replicates', 100);
res.(monk).pro.center   = c;
res.(monk).all.center   = c;
[~, id]                 = sort(cluster_idx);
hist(cluster_idx);
title('cluster distribution');

subplot(3, 3, 2);
imagesc(dat(id,:));
title([monk ' clustered']);
axis off;

% project responses on principal components
[~, loadings] = princomp(dat);
subplot(3, 3, 3)
for j = 1:conf.n_cluster

    id = (cluster_idx == j);
    plot(loadings(id,1), loadings(id,2), ['.' colors{j}]);
    hold on
    xlabel('1st PC');
    ylabel('2nd PC');
end

% kmeans resid test
subplot(3, 3, 4:6);
resid = NaN(1,10);
for j = 1:10
    [~, ~, sumd]  = kmeans(dat, j, 'replicates', 100);
    resid(j)      = sum(sumd);
end
plot(resid)
title('mean field sum')

subplot(3, 3, 7:9);
dat   = dat([1:6 8:end],:);
resid = NaN(1,10);
for j = 1:10
    [~, ~, sumd]  = kmeans(dat, j, 'replicates', 100);
    resid(j)      = sum(sumd);
end
plot(resid)
title('mean field sum (outlier removed)')


saveas(h, [conf.cur_res_fold  'raw_clusters_' monk '.' conf.image_format]);
close(h);
clear c cluster_idx dat h id j loadings monk resid sumd


%% %%%% VEGA
monk    = 'vega';
modi    = {'pro', 'sup', 'all'};
conf.n_cluster = 3;
h       = figure('Visible','off');
for i = 1:length(modi)

    dat                         = res.(monk).(modi{i}).flat;
    [cluster_idx c]             = kmeans(dat, conf.n_cluster, 'replicates', 100);
    res.(monk).(modi{i}).center = c;
    res.(monk).(modi{i}).idx    = cluster_idx;
end

subplot(3, 2, 1)
hist(res.(monk).all.idx);
title('cluster distribution');

subplot(3, 2, 2)
[~, id]          = sort(res.(monk).all.idx);
imagesc(dat(id,:));
title([monk ' clustered ' modi{i}]);
axis off;

% kmeans resid test
subplot(3, 2, 3:4)
resid = NaN(1,10);
for j = 1:10
    [~, ~, sumd]  = kmeans(dat, j, 'replicates', 100);
    resid(j)      = sum(sumd);
end
plot(resid)
title('mean field sum')


% project responses on principal components
[~, loadings] = princomp(res.(monk).all.flat);

subplot(3, 2, 5)
for j = 1:conf.n_cluster

    id = (cluster_idx == j);

    plot(loadings(id,1), loadings(id,2), ['.' colors{j}]);
    hold on
    xlabel('1st PC');
    ylabel('2nd PC');
    title(modi{i});
end

subplot(3, 2, 6)
[a, b, score]   = match_syns(res.(monk).pro.center, res.(monk).sup.center, 1);
tmp(1:2:conf.n_cluster*2-1,:)  = a;
tmp(2:2:conf.n_cluster*2,:)    = b;
imagesc(tmp);
title(num2str(score));

saveas(h, [conf.cur_res_fold  'raw_clusters_' monk '.' conf.image_format]);
close(h);
clear a b c cluster_idx colors dat h i id j loadings modi monk resid score sumd tmp




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
        for i = find(idx.(monk) & idx.(mod{k}))
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
