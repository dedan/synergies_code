

%%
path = '/Volumes/LAB/';
monk = 'vega';

srcdir = [path monk filesep];
load([path 'results' filesep 'data' filesep 'evoked_data_' monk '.mat']);
load([path 'results' filesep 'data' filesep 'nat_mov_res_' monk '.mat']);
load([path 'results' filesep 'data' filesep 'channels.mat']);


%% flexor extensor stuff
% when we have a significant response, is it mostly flexors or extensors or
% doesn't it make a difference

t = [];
for i = 1:length(resps)
    t = [t channels.(monk).type(resps(i).fields)]; %#ok<AGROW>
end

figure(2)
subplot 211
hist(channels.(monk).type)
xlabel(channels.description);
title('muscle types');

subplot 212
hist(t)
title('muscle types of response fields');





% investigation of the pd distribution and of direction of a movement of a
% response according to the pds
n = 6;
figure
subplot(n,1,1:3);
imagesc(normr(vertcat(resps.response)));

subplot(n,1,4);
[x y] = pol2cart(nat_mov_res.pds, ones(1,16));
feather(x,y);

% look only at sessions with a certai fieldsize
interesting = find([resps.field] > 8);

subplot(n,1,5);
r = 0.01 * (resps(interesting(1)).response - min(resps(interesting(2)).response));
% r = normr(resps(interesting(1)).response);

[x y] = pol2cart(nat_mov_res.pds, r);
feather(x,y);

subplot(n,1,6);
feather(sum(x), sum(y));



