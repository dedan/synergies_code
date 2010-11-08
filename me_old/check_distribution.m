
% look at the distribution of the values around stimulations

addpath(genpath('..')); %#ok<MCAP>

window_start   = -6;
window_end     = 30;
sampling_rate  = 500;
stim_value     = 150;

apply_filter   = 0;
rand_stim      = 0;


% folder = '/Volumes/new_disk/darma/';
 folder = '/Volumes/new_disk/vega/';

%% collect information
% get all the responses for a monkey
[res_matrix res_info] = responses(window_start, window_end, sampling_rate, apply_filter, folder, stim_value, rand_stim);
rand_stim             = 1;
[res_matrix_rnd res_info_rnd] = responses(window_start, window_end, sampling_rate, apply_filter, folder, stim_value, rand_stim);

% TODO, das hier normen ?

%% plot the histograms
figure(1);
subplot(2,1,1);
hist(res_matrix(:),100);
title('real stimulation');

subplot(2,1,2);
hist(res_matrix_rnd(:),100);
title('random stimulation');


%% compute means
for i = 1:size(res_matrix)
   norm_res(i,:) = res_matrix(i,:) ./ norm(res_matrix(i,:)); %#ok<AGROW>
   norm_res_rnd(i,:) = res_matrix_rnd(i,:) ./ norm(res_matrix_rnd(i,:)); %#ok<AGROW>
end
norm_res(isnan(norm_res)) = 0;
norm_res_rnd(isnan(norm_res_rnd)) = 0;


%% plot normalized histograms

figure(2);
subplot(2,1,1);
hist(norm_res(:),100);
title('real stimulation norm');

subplot(2,1,2);
hist(norm_res_rnd(:),100);
title('random stimulation norm');


