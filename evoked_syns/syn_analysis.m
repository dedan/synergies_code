
%% add path
addpath(genpath('../lib')); %#ok<MCAP>
addpath(genpath('../preprocess')); %#ok<MCAP>


%% set configuration

clear

config = struct;

config.comment = 'will mehr session info in den fin_res struct schreiben';

config.monk                     = 'chalva';
config.modi                     = {'pro', 'sup'};

% general settings
config.sort_out_zerofield       = 1;                                          % sort out sessions with 0 muscle response field
config.skip_window_computation  = 1;
config.do_resid_test            = 0;
config.debug_images             = 0;
config.inflag                   = true;                                          % write info messages
config.erflag                   = true;                                          % write error messages
config.image_format             = 'pdf';                                     % format of images


% emg channels
config.channels            = [11 12 13 14 21 22 23 24 31 32 33 34 41 42 43 44];             % emg channel names

% response
config.stim_value          = 150;                                          % look only at stimulations around this value

% matrix factorization
config.dimensions          = 3;                                                             % reduce to this number of synergies
config.Niter_res_test      = 10;                                           % number of iterations for the residual test
config.Niter_exploration   = 50;                                          % number of iterations for nmf exploration
config.n_best              = 20;
config.opt                 = statset('MaxIter',5);


% files and folders
config.result_folder       = '~/Documents/uni/yifat_lab/results/';         % path for output


% create folder for output of this session
% config.start_time             = datestr(now, 'ddmmyy_HHMM');
% config.cur_res_fold           = [config.result_folder config.start_time filesep];
% mkdir(config.cur_res_fold);
config.cur_res_fold           = [config.result_folder 'syn_an_debug' filesep];
mkdir(config.cur_res_fold)


% sort sessions again according to the handsorted flags (see also ../response/sort_files)
config.handsorted          = '';

% write all the output into a file
diary([config.cur_res_fold 'info.txt']);
disp(config);
diary('off');
diary([config.cur_res_fold 'out.txt']);

disp(['doing analysis for: ' config.monk]);



%% get responses

if config.skip_window_computation
    disp('skip window calculation, take data from previous session..');
else
    disp('calculate responses .. ');
    runscript4evoked('/Volumes/LAB/', {config.monk})
end
load([config.result_folder 'data' filesep 'evoked_data_' config.monk]);



%% separate and filter responses


disp(' ');
disp(['found ' int2str(length(resps)) ' subsessions in which stimulation took place']);

% sort out the sessions which contain artefacts
load([config.result_folder 'data' filesep config.handsorted 'sort_' config.monk]);
resps = resps(flags ~= 2); 
disp(['sorted out ' int2str(length(find(flags == 2))) ' subsessions because of artefacts (handsorted)']);
disp([int2str(length(resps)) ' subsessions remaining']);


% plot the field size histogram
h = figure('Visible','off');
hist([resps.field],12);
saveas(h, [config.cur_res_fold  'resp_field.' config.image_format]);
close(h);

% sort out the sessions with fieldsize 0
disp(' ');
disp(['sorted out ' int2str(length(find([resps.field] == 0))) ' subsessions because of fieldsize was 0']);
resps = resps([resps.field] ~= 0); 
disp([int2str(length(resps)) ' subsessions remaining']);

% check the stimulation amps for the remaining files
amps = NaN(1,length(resps));
for i = 1:length(resps)
    amps(i) = abs(resps(i).info.amp);
end
% plot the stim amp distribution
h = figure('Visible','off');
hist(amps);
saveas(h, [config.cur_res_fold  'stim_amp_dist.' config.image_format]);
close(h);


% sort the resulting vectors according to the hand position of the subsession
fin_res = sep_results_handposition(resps);

for i = 1:length(fin_res)
    
    disp(' ');
    disp(config.modi{i});
    disp(['number of recordings: ' int2str(length(fin_res(i).dat))]);
end

disp('calculated and separated responses');
clear data filtered_data resp flags sep_results 






%% distribution of final data 

% plot distribution
h = figure('Visible','off');
for i = 1:length(fin_res)
    subplot(length(fin_res),1,i)
    hist(fin_res(i).dat(:),100);
    title(['handposition: ' int2str(i)]);
end
saveas(h, [config.cur_res_fold 'response_dist' '.' config.image_format]);
close(h);





%% resid tests to see to which dimension data can be reduced

if config.do_resid_test
    
    disp(' ');
    disp('doing resid test ..');
        
    % doing the residual test for nmf
    for i = 1:length(fin_res)
        resid_res(i,:)  = test_resid_nmf(fin_res(i).dat, config); 
        
        %repeat for shuffled data as only one random shuffling might not be representative
        tmp = zeros(config.Niter_res_test, min(size(fin_res(i).dat)));
        for j = 1:config.Niter_res_test
            tmp(j,:)    = test_resid_nmf(shuffle_inc(fin_res(i).dat), config);
        end
        resid_resr(i,:) = mean(tmp); %;
        
    end
    
    % add a line to show 10 %
    resid_res(i+1,:)    = ones(1,size(resid_res,2)) * 10;
    resid_resr(i+1,:)   = ones(1,size(resid_resr,2)) * 10;
    
    % residual test for pcaica
    for i = 1:length(fin_res)
        resid_res = [resid_res; test_resid_pcaica(fin_res(i).dat, config)]; %#ok<AGROW>

        %repeat for shuffled data as only one random shuffling might not be representative
        tmp = zeros(config.Niter_res_test, min(size(fin_res(i).dat)));
        for j = 1:config.Niter_res_test
            tmp(j,:)    = test_resid_pcaica(shuffle_inc(fin_res(i).dat), config);
        end
        resid_resr      = [resid_resr; mean(tmp)]; %#ok<AGROW>
    end
    
    first = ones( (length(fin_res)*2) +1,1)*100; 
    first(length(fin_res)+1) = 10;    
    h = createfigure1([first resid_res]', [first resid_resr]',config);
    saveas(h, [config.cur_res_fold 'resid.' config.image_format]);
    close(h);
    
    disp(' ');
    disp('finished resid test !');
    
else
    disp('resid test skipped');
end
clear resid_res resid_resr



%% nmf stuff

disp(' ');
disp('start the search for synergists..');

for i = 1:length(fin_res)
    
    nmf_res     = nmf_explore(fin_res(i).dat, config);
    nmf_res_r   = nmf_explore(shuffle_inc(fin_res(i).dat), config);
    
    h = figure('Visible','off');
    imagesc(nmf_res.flat);
    title(['standard deviation of group size: ' num2str(nmf_res.std)]);
    saveas(h, [config.cur_res_fold  'nmf_expl_stab' int2str(i) '.' config.image_format]);
    close(h);
    disp(['standard deviation of group size: ' num2str(nmf_res.std)]);
    
    fin_res(i).nmf_syns  = nmf_res.syns;
    fin_res(i).nmf_res   = nmf_res;
    fin_res(i).nmf_res_r = nmf_res_r;
end

disp(' ');
%clear nmf_res



%% syn analysis with pcaica

for i = 1:length(fin_res)
    Wpi                  = pcaica(fin_res(i).dat, config.dimensions);
    fin_res(i).pca_syns  = Wpi';
end;
disp('finished the search for synergists !');
clear Wpi


%% compare the three best nmf synergies with the three best pcaica

for i = 1:length(fin_res)
    
    h = figure('Visible','off');
    nmf_sorting = syn_compare(fin_res(i).pca_syns, fin_res(i).nmf_syns, 'pca', 'nmf');
    saveas(h, [config.cur_res_fold  'best_syn_relation' int2str(i) '.' config.image_format]);
    close(h);
    
    % sort the nmf synergies according to the pca synergies
    fin_res(i).nmf_syns = fin_res(i).nmf_syns(nmf_sorting,:);
    
end
clear nmf_sorting


%% compare the three nmf synergists for pronation and supination

if length(fin_res) > 1
    h = figure('Visible','off');
    syn_compare(fin_res(1).nmf_syns, fin_res(2).nmf_syns, 'pronation', 'supination');
    saveas(h, [config.cur_res_fold  'between_pos_syn_relation.' config.image_format]);
    close(h);
    
    h = figure('Visible','off');
    syn_compare(fin_res(1).pca_syns, fin_res(2).pca_syns, 'pronation', 'supination');
    saveas(h, [config.cur_res_fold  'between_pos_syn_relation_pca.' config.image_format]);
    close(h);
end

%% save configuration
save([config.cur_res_fold 'fin_res'], 'fin_res');
save([config.cur_res_fold 'config'], 'config' );
diary('off');

disp('finished !');
