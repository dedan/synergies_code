
clear

path  = '/Volumes/LAB/results/data/';

conf = struct;
conf.Niter_exploration = 50;
conf.n_best            = 20;
conf.opt               = statset('MaxIter',5);
conf.dimensions        = 3;


% this will load a struct called sessions
load([path 'all_data_vega']);
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



%% do we find the correct synergies also when the order of the model is
% larger than it should be or has to be?
figure(3)
config = struct;
config.Niter_exploration= 30;
config.n_best           = 30;
config.opt              = statset('MaxIter',50);

dat = sessions(1).mats(1).data_raw;

config.dim  = 3;
res3        = nmf_explore(dat, config);

config.dim  = 4;
res4        = nmf_explore(dat, config);

[m1, m2, scores] = match_syns(res3.syns, res4.syns, 1);
subplot 211; imagesc(m1); subplot 212; imagesc(m2)
title(['scores: ' num2str(scores)]);



%% gleiche ergebnisse ob normalisiert oder nicht?
figure(4)
config = struct;
config.Niter_exploration= 30;
config.n_best           = 30;
config.opt              = statset('MaxIter',50);
config.dim  = 3;

res        = nmf_explore(sessions(1).mats(1).data_raw, config);
res_norm   = nmf_explore(normr(sessions(1).mats(1).data_raw), config);

[m1, m2, scores] = match_syns(res.syns, res_norm.syns, 1);
subplot 211; imagesc(m1); subplot 212; imagesc(m2)
title(['scores: ' num2str(scores)]);