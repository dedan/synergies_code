% in this file i wanted to compare the synergies from evoked responses with
% the synergies from natural movement. Finally I did it with the projection
% as I did not find stable synergies for evoked data for a long time and
% the projection does not rely on this. Now, as I also found stable
% synergies for the evoked conditions it should also be possible to compare
% the subspaces by looking at the principal angles between the spaces.

% so this script is not in a useful version, but could be used as a start..

%% load data

res_folder  = '~/Documents/uni/yifat_lab/results/';
load([res_folder 'nonevoked_syns/nonevoked_results.mat']);
load('~/Documents/uni/yifat_lab/results/data/060909_2343.mat');

e_base = struct;
e_base(1).s = fin_res(1).nmf_syns';
e_base(2).s = fin_res(2).nmf_syns';
e_base(3).s = fin_res(3).nmf_syns';

n_base = struct;
n_base(1).s = nonevoked_res.fin_syns.nmf_all_pro';
n_base(2).s = nonevoked_res.fin_syns.nmf_all_sup';
n_base(3).s = nonevoked_res.fin_syns.nmf_all';

e_base(3).s
subspace(n_base(3).s,n_base(3).s + randn(size(n_base(3).s)))



%% compute baseline for normalization

% the normalization I used here is from the M. Tresch 2006 Paper about
% different matrix factorization techniques
val = NaN(3,10000);
for i = 1:10000
    for j = 1:3
        perms = randperm(length(fin_res(j).dat));
        val(j,i) = subspace(fin_res(j).dat(perms(1:3),:)', nonevoked_res.fin_syns.nmf');
    end 
end

figure(1)
for i = 1:3
    subplot(3,1,i)
    hist(val(i,:),100);
    xlim([0 1.6]);
end

d_b = mean(val,2);

for i = 1:3
    d = subspace(e_base(i).s, n_base(i).s);
    d_n = (d-d_b(i))/(0 - d_b(i));
    disp(d_n);
end








