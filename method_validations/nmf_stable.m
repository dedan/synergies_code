
clear

path  = '~/Documents/uni/yifat_lab/synergies_code/test_data/';

conf = struct;
conf.Niter_exploration = 50;
conf.n_best            = 20;
conf.opt               = statset('MaxIter',5);
conf.dimensions        = 3;


% this will load a struct called sessions
load([path 'all_data_syn']);
addpath('../lib'); 


%% susceptibility to outliers
% how much does an outlier in the data_raw matrix influence the outcome of
% nmf

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



%% stability of nmf method
% because it is a validation of the method and not really data analysis
% this investigation is moved to this separate script. 
% quantification of nmf stability is done by looking at the std of indeces
% which assign synergies from the several runs to protoypes of a kmeans
% clustering of the synergies from several runs

figure(2)
hist([sessions.syn_pro_std]);

