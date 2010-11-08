
%% add path
addpath(genpath('..')); %#ok<MCAP>

%% set configuration

config = struct;
True   = 1;
False  = 0;

config.comment = 'kmeans oft, einen weg sortierung';

config.skip_window_computation   = 0;
config.do_resid_test             = 1;

config.explore_and_cluster       = 1;

% debugging stuff
config.n_debug_images         = 10;                                        % number of debug images
config.debug_sessions         = [];                                        % variable to be filled with the sessions that are used for debug
config.debug_images           = False;                                      % print debug images
config.inflag                 = False;                                     % write info messages
config.erflag                 = True;                                      % write error messages
config.image_format           = 'tiff';                                    % format of error messages
config.result_folder          = '~/Documents/uni/yifat_lab/results/';      % path for output

% create folder for output of this session
d = clock;
config.start_time             = [int2str(d(3)) int2str(d(2)) int2str(d(1)) '__' int2str(d(4)) int2str(d(5))];
config.cur_res_fold           = [config.result_folder config.start_time filesep];
mkdir(config.cur_res_fold);


% 1: without, 2: std, 3: von mises, 4: von mises++, 5: maximum, 6: not emg
% only in window from baseline
config.normalization       = 6;
config.dimensions          = 10;                                           % dimensions for synergist analysis
config.channels2take       = logical([0  1  1  1  1  1  1  1  1  1  1  1  0  0  0  0 ]); %channel 1 and 13 to 16 are ignored
config.channels            = [11 12 13 14 21 22 23 24 31 32 33 34 41 42 43 44];     % emg channel names
config.n_channels          = length(config.channels);                      % number of emg channels
config.n_used_channels     = length(find(config.channels2take));
config.window              = [-30 30];                                     % start and end of average window
config.int_window          = [6 20];                                       % integrate window only in this part
config.sampling_rate       = 10000;                                          % downsampling to this rate (acts as lowpass)
config.apply_filter        = False;                                         % artefact filter around stimulation time
config.stim_value          = 150;                                          % look only at stimulations around this value
config.Niter_res_test      = 50;                                           % number of iterations for the residual test
config.Niter_exploration   = 200;                                          % number of iterations for nmf exploration
config.rand_stim           = False;                                        % get average windows around random stimulation times
config.monkey              = 'vega';                                       % monkey to analyze
config.volume              = '/Volumes/lab/';                              % harddisk which contains the data
config.folder              = [config.volume config.monkey filesep];
config.pd_folder           = [config.folder 'pd_files' filesep];           % PD files (preferred direction)
config.dat_folder          = [config.folder 'data' filesep];               % folder with data



diary([config.cur_res_fold 'info.txt']);
disp(config);
diary('off');
diary([config.cur_res_fold 'out.txt']);


% config consistency
if config.normalization == 6 && config.window(1) ~= -30
   error('this normalization needs the prestimulus part of the signal');

end


%% get responses and normalize

% struct for the final results
fin_res = struct;

disp('calculate responses and normalize.. ');

if config.skip_window_computation
   disp('skip window calculation, take data from previous session..');
else

   % get all subessions in which a stimulation took place
   data = get_all_stimulations(config.dat_folder);

   % filter subessions according to StimAmp
   filtered_data = stimulations_at(data, config.stim_value);

   % set which sessions are taken for debug images
   perms                   = randperm(length(filtered_data));
   config.debug_sessions   = [filtered_data(perms(1:config.n_debug_images)).id];
   disp(config.debug_sessions);

   % get all the responses for a monkey
   [res_matrix res_info] = responses(filtered_data, config);

   % normalize the vector length for each session over all channels to 1
   norm_res = norm_rows(res_matrix(:,config.channels2take));

   disp('calculated responses and normalized them !');
end


%% separation of results

% sort the resulting vectors according to the hand position of the subsession
sep_results = sep_results_handposition(res_info, norm_res);


% NOTE because for vega the mapping of the muscles changed after session 25
% the results have to be separated again
remap_idx = struct;
for i = 1:length(sep_results)
   if strcmp(config.monkey, 'vega')
      remap_idx(i).before  = find([sep_results(i).info.id] < 26);
      remap_idx(i).after   = find([sep_results(i).info.id] > 25);
   else
      remap_idx(i).before  = [];
      remap_idx(i).after   = [sep_results(i).info.id];
   end
end

disp('separated responses');

% NOTE analysis from here only for after remapping because it exists much more
% data and also the PDs are computed only for after remapping
for i = 1:length(sep_results)
   fin_res(i).dat = sep_results(i).dat(remap_idx(i).after,:);
   fin_res(i).info = sep_results(i).info(remap_idx(i).after);
end



%% analyze variance in groups

if config.do_resid_test

   disp(' ');
   disp('doing resid test ..');

   % doing the residual test for nmf
   for i = 1:length(fin_res)
      resid_res(i,:) = [100 test_resid_old(fin_res(i).dat, config.Niter_res_test)]; %#ok<AGROW>
   end

   % add a line to show 10 %
   resid_res(i+1,:) = ones(1,size(resid_res,2)) * 10;

   % doing the residual test for pcaica
   for i = 1:length(fin_res)
      resid_res_pcaica(i,:) = [100 test_resid_pcaica(fin_res(i).dat, config)]; %#ok<AGROW>
   end

   % add it to othe resid test data and plot
   resid_res = [resid_res; [resid_res_pcaica zeros(2,size(resid_res,2)-size(resid_res_pcaica,2))]];
   h = createfigure(resid_res');
   saveas(h, [config.cur_res_fold 'resid.' config.image_format]);
   close(h);

   disp(' ');
   disp('finished resid test !');

else
   disp('resid test skipped');
end


%% exploration nmf space, search for global maximum

disp(' ');
disp('start the search for synergists..');

for j = 1:length(fin_res)

   if config.explore_and_cluster

      dim            = config.dimensions;
      iter           = config.Niter_exploration;

      fin_res(j).all_patterns   = zeros(dim*iter,config.n_used_channels);

      % exploration
      for i = 1:iter
         W = NMF(fin_res(j).dat',dim)';
         fin_res(j).all_patterns(dim*i - dim +1:dim*i,:) = W;
      end

      % find prototypes by k_means clustering
      min = 100000000000;
      for i = 1:config.Niter_exploration
         [IDX,protos, err]      = kmeans(fin_res(j).all_patterns,dim);
         if sum(err) < min
            min = sum(err);
            fin_res(j).protos   = protos;
            fin_res(j).idx      = IDX;
         end
      end

      % get the loadings for the found prototypes
      [W H score] = NMF(fin_res(j).dat', config.dimensions, fin_res(j).protos');
      fin_res(j).loadings = H;

      h = figure('Visible','off');
      hist(fin_res(j).idx, config.dimensions);
      saveas(h, [config.cur_res_fold  'kmeans_dist' int2str(j) '.' config.image_format]);
      close(h);

   else
      % simple exploration without clustering
      min = 100000000000;
      for i = 1:config.Niter_exploration
         [W,H, score] = NMF(fin_res(j).dat', config.dimensions);
         if score < min
            min = score;
            fin_res(j).protos = W';
            fin_res(j).loadings = H;
         end
      end
   end
end

disp(' ');
disp('finished the search for synergists !');


%% look for the best synergists (sort synergists)
for i = 1:length(fin_res)

   res = ones(config.dimensions,1);
   for j = 1:config.dimensions
      choose      = false(config.dimensions,1);
      choose(j)   = true;
      res(j)      = sum(sum(abs(fin_res(i).dat' - fin_res(i).protos(choose,:)' * fin_res(i).loadings(choose,:))));
   end

   [sorted idx]         = sort(res, 'descend');

   fin_res(i).errors    = sorted;
   fin_res(i).syn_order = idx;
end


%% syn analysis with pcaica

for i = 1:length(fin_res)
   [Wpi,loadings, score]   = pcaica(fin_res(i).dat, config.dimensions);
   fin_res(i).pca_syns     = Wpi';
end;




%% compute prefered directions for synergists


disp(' ');
disp('computing preferred directions for synergists..');

pds = get_pds(config);

fin_res(1).channel_pd = pds(1,:);
fin_res(2).channel_pd = pds(2,:);

disp(' ');
disp('computed pds for synergists !');


%% synergy comparison plots

for i = 1:length(fin_res)
   
   % norm the prototypes and sort the nmf protos
   normed_protos = norm_rows(fin_res(i).protos);
   normed_protos = normed_protos(fin_res(i).syn_order,:);
   normed_pcaica = norm_rows(fin_res(i).pca_syns);
   
   % map the most similar protos from both procedures
   mapping              = map(normed_protos, normed_pcaica);
   fin_res(i).mapping   = mapping;

   n_best = 5;
   h = figure('Visible','off');
   for j = 1:n_best
      subplot(4,n_best,j);
      bar(normed_protos(j,:));
   end
   for j = 1:n_best
      subplot(4,n_best,n_best+j);
      bar(normed_pcaica(j,:));
   end
   for j = 1:n_best
      subplot(4,n_best,2*n_best+j);
      bar(normed_protos(mapping(j,2),:));
   end
   for j = 1:n_best
      subplot(4,n_best,3*n_best+j);
      bar(normed_pcaica(mapping(j,3),:));
   end
   saveas(h, [config.cur_res_fold 'syn_comparison' int2str(i) '.' config.image_format]);
   close(h);

   h = figure('Visible','off');
   hold on;
   for j = 1:3
      plot(normed_protos(mapping(j,2),:), normed_pcaica(mapping(j,3),:), '.');
   end
   hold off;
   saveas(h, [config.cur_res_fold 'syn_value_dist' int2str(i) '.' config.image_format]);
   close(h);
end




%% create final plots


disp(' ');
disp('plotting the final results..');

for j = 1:length(fin_res)

   h = figure('Visible','off');
   subplot(4,1,1);
   imagesc(fin_res(j).protos(fin_res(j).syn_order,:));
   title('protos ordered');

   subplot(4,1,2);
   imagesc(fin_res(j).pca_syns);
   title('pcaica synergists');

   subplot(4,1,3:4);
   plot(fin_res(j).errors);
   saveas(h, [config.cur_res_fold 'final' int2str(j) '.' config.image_format]);
   close(h);

   h = figure('Visible','off');
   % rose plots
   sorted_syns = fin_res(j).protos(fin_res(j).syn_order,:);
   for k = 1:4
      subplot(2,2,k);
      rose_agg = [];
      used_pds = fin_res(j).channel_pd(config.channels2take);
      for i = 1:config.n_used_channels
         rose_agg = [rose_agg ones(1,floor(sorted_syns(k,i) * 100)) * used_pds(i)]; %#ok<AGROW>
      end
      h_fake = rose(ones(1,500));
      hold on;
      rose(rose_agg);
      title(int2str(fin_res(j).syn_order(k)));
      set(h_fake, 'Visible', 'Off');
   end
   saveas(h, [config.cur_res_fold 'final_rose' int2str(j) '.' config.image_format]);
   close(h);

end

% plot distribution
h = figure('Visible','off');
both_hand = [fin_res(1).dat(:); fin_res(2).dat(:)];
both_hand = both_hand(both_hand ~= 0);
hist(both_hand, 100);
saveas(h, [config.cur_res_fold 'distribution' '.' config.image_format]);
close(h);

disp(' ');
disp('plotted the final results !');




%% report

disp(' ');
disp('Summary: ');

disp(' ');
disp('Pronation: ');
disp(['number of recordings: ' int2str(length(fin_res(1).dat))]);

disp(' ');
disp('Supination: ');
disp(['number of recordings: ' int2str(length(fin_res(2).dat))]);


%% save configuration
fin_res(1).config = config;
fin_res(2).config = config;
save([config.result_folder 'data' filesep config.start_time], 'fin_res');
save([config.cur_res_fold 'data'], 'fin_res' );
diary('off');

disp('finished !');
