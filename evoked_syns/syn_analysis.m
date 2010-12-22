
%% add path
addpath(genpath('../lib')); %#ok<MCAP>
addpath(genpath('../preprocess')); %#ok<MCAP>


%% set confuration

clear

conf = struct;

conf.comment = 'will mehr session info in den fin_res struct schreiben';

conf.monks                    = {'chalva'};
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
conf.dimensions          = 3;             % reduce to this number of synergies
conf.Niter_res_test      = 10;            % number of iterations for the residual test
conf.Niter_exploration   = 50;            % number of iterations for nmf exploration
conf.n_best              = 20;            % take only n_best explorations
conf.opt                 = statset('MaxIter',50);     % number of early explorations


% files and folders
conf.result_folder       = '~/Documents/uni/yifat_lab/results/';
conf.inpath              = '/Volumes/LAB/results/';
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
clear tmp;


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
    load([conf.result_folder 'data' filesep conf.handsorted 'sort_' monk]);
    
    % wow, what a line. I do this weird indexing to make it the same length
    % as flags
    idx.(monk)(idx.(monk)) = idx.(monk)(idx.(monk)) & (flags ~= 2);
    disp(['sorted out ' int2str(length(find(flags == 2))) ...
        ' subsessions because of artefacts (handsorted)']);
    disp([int2str(length(find(idx.(monk)))) ' subsessions remaining']);


    % plot the field size histogram
    figure(h1)
    subplot(length(conf.monks), 1, m)
    hist([resps(idx.(monk)).field],0:length(conf.channels));
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
    
    res.(monk).pro.flat = vertcat(resps(idx.(monk) & idx.pro).response);
    res.(monk).sup.flat = vertcat(resps(idx.(monk) & idx.sup).response);
    if ~isempty(res.(monk).sup.flat)
        res.(monk).all.flat = vertcat(res.(monk).pro.flat, res.(monk).sup.flat);
    else
        res.(monk).all.flat = [];
    end
    
    disp(['pronation - number of recordings: ' num2str(length(find(idx.(monk) & idx.pro)))]);
    disp(['supination - number of recordings: ' num2str(length(find(idx.(monk) & idx.sup)))]);
    
end

saveas(h1, [conf.cur_res_fold  'resp_fields.' conf.image_format]);
close(h1);
saveas(h2, [conf.cur_res_fold  'stim_amp_dist.' conf.image_format]);
close(h2);

disp(' ');
disp('calculated and separated responses');
clear flags monk h1 h2





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





%% distribution of final data 

% plot distribution
h = figure('Visible','off');

for m = 1:length(conf.monks)
    monk = char(conf.monks{m});
    
    subplot(2, length(conf.monks), m)
    [n,xout] = hist([resps(idx.(monk) & idx.pro).response], 100);
    n        = (n / max(n)) *100;
    bar(xout,n)
    title([monk ' pronation']);
    
    subplot(2, length(conf.monks), m + length(conf.monks))
    [n,xout] = hist([resps(idx.(monk) & idx.sup).response], 100);
    n        = (n / max(n)) *100;
    bar(xout,n)
    title([monk ' supination']);
end    
saveas(h, [conf.cur_res_fold 'response_dist' '.' conf.image_format]);
close(h);




%% resid tests to see to which dimension data can be reduced

if conf.do_resid_test
    
    disp(' ');
    disp('doing resid test ..');
    
    for m = 1:length(conf.monks)
        
        monk = char(conf.monks{m});
        clear resid_res;
                
        
        tmp_pro             = res.(monk).pro.flat;
        resid_res(1,:)      = test_resid_nmf(tmp_pro, conf); 

        %repeat for shuffled data as only one random shuffling might not be representative
        tmp = zeros(conf.Niter_res_test, min(size(tmp_pro)));
        for j = 1:conf.Niter_res_test
            tmp(j,:)    = test_resid_nmf(shuffle_inc(tmp_pro), conf);
        end
        resid_res(2,:) = mean(tmp); 
        l = {'resid pro', 'shuffled resid pro'};

        % if supination data available
        tmp_sup             = res.(monk).sup.flat;
        if ~isempty(tmp_sup)
            resid_res(3,:)  = test_resid_nmf(tmp_sup, conf); 
            tmp = zeros(conf.Niter_res_test, min(size(tmp_sup)));        
            for j = 1:conf.Niter_res_test
                tmp(j,:)    = test_resid_nmf(shuffle_inc(tmp_sup), conf);
            end
            resid_res(4,:)  = mean(tmp);
            
            tmp_all         = res.(monk).all.flat;
            resid_res(5,:)  = test_resid_nmf(tmp_all, conf);
            tmp = zeros(conf.Niter_res_test, min(size(tmp_all)));
            for j = 1:conf.Niter_res_test
                tmp(j,:)    = test_resid_nmf(shuffle_inc(tmp_all), conf);
            end
            resid_res(6,:) = mean(tmp);

            l = {'resid pro', 'shuffled resid pro', 'resid sup', ...
                'shuffled resid sup', 'resid all', 'shuffled resid all'};
        end
        
        % add a line to show 10 %
        resid_res(size(resid_res,1)+1,:)    = ones(1,size(resid_res,2)) * 10; %#ok<SAGROW>
        first                               = ones(size(resid_res,1),1)*100;
        first(size(resid_res,1))            = 10;
        resid_res                           = [first resid_res]; %#ok<AGROW>
        
        h = figure('Visible','off');
        plot(0:size(resid_res,2)-1, resid_res');
        legend(l);
        saveas(h, [conf.cur_res_fold 'resid_' monk '.' conf.image_format]);
        close(h);        
    end
    
    disp(' ');
    disp('finished resid test !');
    
else
    disp('resid test skipped');
end
clear resid_res tmp_pro tmp_sup tmpp tmps first



%% matrix factorizations

disp(' ');
disp('start the search for synergists..');

pos = {'pro', 'sup', 'all'};

for m = 1:length(conf.monks)    
    monk = char(conf.monks{m});
    h = figure('Visible','off');
    
    for mo = 1:length(pos)
        mod = char(pos{mo});
        
        if isempty(res.(monk).(mod).flat)
            res.(monk).(mod).nmf = [];
            res.(monk).(mod).pca = [];
            continue;
        end
        % compute synergies
        nmf_res      = nmf_explore(res.(monk).(mod).flat, conf);
        disp(['standard deviation of group size: ' num2str(nmf_res.std)]);
        pca_res      = pcaica(res.(monk).(mod).flat, conf.dim)';
        
        % norm, match and sort them
        [pca_res, nmf_res.syns, scores] = match_syns(pca_res, nmf_res.syns);
        
        % plot
        for i = 1:conf.dim
            subplot(length(pos), conf.dim+1, (conf.dim+1)*(mo-1) + i)
            bar( [pca_res(i,:)' nmf_res.syns(i,:)']);
            axis off
            title(['#' int2str(i) ' sc: ' num2str(scores(i))]);
        end
        subplot(length(pos), conf.dim+1, conf.dim*mo +1)
        plot(pca_res(:), nmf_res.syns(:), '.');
        [p, r] = corrcoef(pca_res(:), nmf_res.syns(:));
        title([num2str(p(2,1)) ' -- ' num2str(r(2,1))]);

        % save results
        res.(monk).(mod).nmf = nmf_res.syns;    
        res.(monk).(mod).pca = pca_res;    
    end
    saveas(h, [conf.cur_res_fold  'synergies_' monk '.' conf.image_format]);
    close(h);
end

disp(' ');
clear nmf_res pca_res monk mod scores ind



%% save configuration
evoked_res = res;
save([conf.inpath 'data' filesep 'evoked_res'], 'evoked_res');

disp('finished !');
