
clear

outpath  = '~/Documents/uni/yifat_lab/results/natural_mov/';

conf = struct;
conf.Niter_exploration = 50;
conf.opt               = statset('MaxIter',5);
conf.dimensions        = 5;

load([outpath 'all_data_dummy']);
addpath('../lib'); 


%% susceptibility to outliers

n_largest = 10;

% this is a session containing an outlier
i_outl = find([sessions.id] == 46);

dat_out     = sessions(i_outl).mats(1).data_raw;
dat_without = sessions(i_outl).mats(1).data_raw;

sorted = sort(dat_out(:), 'descend');
maxes  = sorted(1:n_largest);

for i = 1:n_largest
    dat_without(dat_out == maxes(i)) = mean(dat_out(:));
end

figure(1)
subplot(5,1,1:3)
imagesc(dat_out)

all_vega = strmatch(sessions(i_outl).monk, {sessions.monk});
all_chan = vertcat(sessions(all_vega).channels);
c2take   = all(all_chan);
                    
nmf_res = nmf_explore(dat_out(:,c2take), conf);

subplot 514
imagesc(nmf_res.syns);

nmf_res = nmf_explore(dat_without(:,c2take), conf);

subplot 515
imagesc(nmf_res.syns);


% if config.nmf_stab_plots == 1
%             h = figure('Visible','off');
%             imagesc(nmf_res.flat);
%             title(['standard deviation of group size: ' num2str(nmf_res.std)]);
%             saveas(h, [config.outpath  'nmf_expl_stab' int2str(i) '_' int2str(j) '.' config.image_format]);
%             close(h);
%         end
% 
%     % show stability
%     disp(['standard deviation of group size: ' num2str(nmf_res.std)]);
%     
%     
%     can be tested with for example 
% (2)session id: 46 - monk: vega
% 
% how much does an outlier in the data_raw matrix influence the outcome of nmf or pca


