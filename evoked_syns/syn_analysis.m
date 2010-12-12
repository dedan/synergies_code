
%% add path
addpath(genpath('../lib')); %#ok<MCAP>
addpath(genpath('../response')); %#ok<MCAP>


%% set configuration

config = struct;

config.comment = 'will mehr session info in den fin_res struct schreiben';

mymap = [linspace(145/255,178/255,32)' linspace(167/255,213/255,32)' linspace(216/255,111/255,32)'];
mymap = [mymap; linspace(178/255,251/255,32)' linspace(213/255,147/255,32)' linspace(111/255,24/255,32)'];
config.mymap = mymap;

% general settings
config.normalize_channels       = 0;                                          % normalize activation within channels
config.sort_out_zerofield       = 1;                                          % sort out sessions with 0 muscle response field
config.skip_window_computation  = 1;
config.do_resid_test            = 0;
config.debug_images             = 0;
config.inflag                   = 0;                                          % write info messages
config.erflag                   = 1;                                          % write error messages
config.image_format             = 'pdf';                                     % format of images


% emg channels
config.channels2take       = logical([0  1  1  1  1  1  1  1  1  1  1  1  0  0  0  0 ]);    %channel 1 and 13 to 16 are ignored
config.channels            = [11 12 13 14 21 22 23 24 31 32 33 34 41 42 43 44];             % emg channel names
config.n_channels          = length(config.channels);                                       % number of emg channels
config.n_used_channels     = length(find(config.channels2take));                            % how many channels are used

% response
% how to normalize the EMG data in the average response window:
% 1: mean activity change, 2: normalization with std
% both normalizations will be computed in response function, which one is
% used is then decided in the post_remapping function
config.normalization       = 1;
config.window              = [-20 20];                                     % start and end of average window
config.int_window          = [6 20];                                       % integrate window only in this part (final response)
config.sampling_rate       = 10000;                                        % downsampling to this rate (acts as lowpass)
config.apply_filter        = 0;                                            % artefact filter around stimulation time
config.stim_value          = 150;                                          % look only at stimulations around this value

% matrix factorization
config.dimensions          = 3;                                                             % reduce to this number of synergies
config.Niter_res_test      = 20;                                           % number of iterations for the residual test
config.Niter_exploration   = 100;                                          % number of iterations for nmf exploration
config.opt                 = statset('MaxIter',5);


% files and folders
config.monkey              = 'vega';                                       % monkey to analyze
config.volume              = '/Volumes/DEDAN_DISK/';                       % harddisk which contains the data
config.folder              = [config.volume config.monkey filesep];
config.pd_folder           = [config.folder 'pd_files' filesep];           % PD files (preferred direction)
config.dat_folder          = [config.folder 'data' filesep];               % folder with data
config.result_folder       = '~/Documents/uni/yifat_lab/results/';         % path for output


% create folder for output of this session
config.start_time             = datestr(now, 'ddmmyy_HHMM');
config.cur_res_fold           = [config.result_folder config.start_time filesep];
mkdir(config.cur_res_fold);


% sort sessions again according to the handsorted flags (see also ../response/sort_files)
config.handsorted          = [config.result_folder 'responses/handsorted_new.mat'];  % handsorted_strict.mat

% write all the output into a file
diary([config.cur_res_fold 'info.txt']);
disp(config);
diary('off');
diary([config.cur_res_fold 'out.txt']);






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
disp(['found ' int2str(length(resp)) ' subsessions in which stimulation took place']);

% sort out the sessions which contain artefacts
load(config.handsorted);
resp = resp(find(flags ~= 2));  %#ok<FNDSB>
disp(['sorted out ' int2str(length(find(flags == 2))) ' subsessions because of artefacts (handsorted)']);
disp([int2str(length(resp)) ' subsessions remaining']);


% plot the field size histogram
h = figure('Visible','off');
hist([resp.field],12);
saveas(h, [config.cur_res_fold  'resp_field.' config.image_format]);
close(h);

% sort out the sessions with fieldsize 0
if config.sort_out_zerofield == 1
    disp(' ');
    disp(['sorted out ' int2str(length(find([resp.field] == 0))) ' subsessions because of fieldsize was 0']);
    resp = resp(find([resp.field] ~= 0)); %#ok<FNDSB>
    disp([int2str(length(resp)) ' subsessions remaining']);
end

% check the stimulation amps for the remaining files
amps = NaN(1,length(resp));
for i = 1:length(resp)
    amps(i) = abs(resp(i).info.amp);
end
% plot the stim amp distribution
h = figure('Visible','off');
hist(amps);
saveas(h, [config.cur_res_fold  'stim_amp_dist.' config.image_format]);
close(h);



% sort the resulting vectors according to the hand position of the subsession
sep_results = sep_results_handposition(resp);

% sort out pre remapping sessions
% NOTE in here is chosen which normalization to take
fin_res = post_remapping(sep_results, config);


disp(' ');
disp('Pronation: ');
disp(['number of recordings: ' int2str(length(fin_res(1).dat))]);

disp(' ');
disp('Supination: ');
disp(['number of recordings: ' int2str(length(fin_res(2).dat))]);


disp('calculated and separated responses');
clear data filtered_data resp flags sep_results 


%% distribution of final data and normalization within channels
% only interesting to chose whether to normalize within channels by max or
% median or mean. However in the ende we do not normalize at all

% plot distribution
h = figure('Visible','off');
for i = 1:length(fin_res)
    subplot(3,1,i)
    hist(fin_res(i).dat(:),100);
    title(['handposition: ' int2str(i)]);
end
saveas(h, [config.cur_res_fold 'distribution' '.' config.image_format]);
close(h);

if config.normalize_channels == 1
    for i = 1:length(fin_res)
        for j = 1:size(fin_res(i).dat,2)
            % NOTE I normalize by mean as max does not seem to be a good
            % solution, when looking at the distribution of the final data.
            fin_res(i).dat(:,j) = fin_res(i).dat(:,j) ./ mean(fin_res(i).dat(:,j));
        end
    end
end




%% resid tests to see to which dimension data can be reduced

if config.do_resid_test
    
    disp(' ');
    disp('doing resid test ..');
        
    % doing the residual test for nmf
    for i = 1:length(fin_res)
        resid_res(i,:)  = test_resid_nmf(fin_res(i).dat, config); %#ok<AGROW>
        
        %repeat for shuffled data as only one random shuffling might not be representative
        tmp = zeros(config.Niter_res_test, min(size(fin_res(i).dat)));
        for j = 1:config.Niter_res_test
            tmp(j,:)    = test_resid_nmf(shuffle_inc(fin_res(i).dat), config);
        end
        resid_resr(i,:) = mean(tmp); %#ok<AGROW>;
        
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
    
    first = ones(7,1)*100; first(4) = 10;    
    h = createfigure1([first resid_res]', [first resid_resr]',config);
    saveas(h, [config.cur_res_fold 'resid.' config.image_format]);
%    close(h);
    
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

        % TODO brauche ich die nmf exploration überhaupt noch, ist das
        % nicht schon sehr stabil?
        nmf_res     = nmf_explore(fin_res(i).dat, config);
        nmf_res_r   = nmf_explore(shuffle_inc(fin_res(i).dat), config);
        
        h = figure('Visible','off');
        imagesc(nmf_res.flat);
        colormap(config.mymap);
        title(['standard deviation of group size: ' num2str(nmf_res.std)]);
        saveas(h, [config.cur_res_fold  'nmf_expl_stab' int2str(i) '.' config.image_format]);
        close(h);
        disp(['standard deviation of group size: ' num2str(nmf_res.std)]);
        
        fin_res(i).nmf_syns = nmf_res.syns;
        fin_res(i).nmf_res  = nmf_res;
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

h = figure('Visible','off');
syn_compare(fin_res(1).nmf_syns, fin_res(2).nmf_syns, 'pronation', 'supination');
saveas(h, [config.cur_res_fold  'between_pos_syn_relation.' config.image_format]);
close(h);

h = figure('Visible','off');
syn_compare(fin_res(1).pca_syns, fin_res(2).pca_syns, 'pronation', 'supination');
saveas(h, [config.cur_res_fold  'between_pos_syn_relation_pca.' config.image_format]);
close(h);


%% save configuration
fin_res(1).config = config;
fin_res(2).config = config;
save([config.result_folder 'data' filesep config.start_time], 'fin_res');
save([config.result_folder 'data' filesep 'latest_evoked'], 'fin_res');
save([config.cur_res_fold 'data'], 'fin_res' );
diary('off');

disp('finished !');
