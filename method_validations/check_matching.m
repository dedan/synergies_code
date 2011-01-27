
% file to test the matching algorithm. Yifat asked me why I match the
% vectors by using the dot product and compute the matching score by
% corrcoef. I don't remember anymore. So I will check the differences in
% this file

clear

%% first create some artifial data

noise_level = 0.2;
syns = rand(3, 11);
ww1 = normr(syns + randn(3, 11) * noise_level);
ww2 = normr(syns + randn(3, 11) * noise_level);


%% or load some test data 
% this data contains some vectors where the dot product gives the correct
% result, the corrcoef not
load('../test_data/dot_corrcoef_matching.mat');


%% plot the original data
figure(1)
subplot 311
imagesc(vertcat(ww1, ww2));


%% use dot product matching
[s1 s2] = match_syns(ww1, ww2, 1);
subplot 312
agg(1:2:5,:) = s1;
agg(2:2:6,:) = s2;
imagesc(agg);


%% use corrcoef matching
[s1 s2] = match_syns(ww1, ww2, 2);
subplot 313
agg(1:2:5,:) = s1;
agg(2:2:6,:) = s2;
imagesc(agg);
