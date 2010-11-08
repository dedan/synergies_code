
% i used this file to create some special plots for the poster..

% my own colormap :-)
mymap = [linspace(145/255,178/255,32)' linspace(167/255,213/255,32)' linspace(216/255,111/255,32)'];
mymap = [mymap; linspace(178/255,251/255,32)' linspace(213/255,147/255,32)' linspace(111/255,24/255,32)'];


%% plot of typical muscle activation in natural movement
load(['~/Documents/uni/yifat_lab/results/nonevoked_syns/' 'all_nonevoked_sessions']);
for i = 1:length(all_sessions)
    imagesc(all_sessions(i).data.mathand(1).data)
    axis image off
    colormap(mymap)
    pause
end

%%
for i = 1:10
    figure
for j = 1:8
    subplot(8,1,j);
    bar(all_sessions(i).data.mathand(1).data(j,:));
end

end

%%
for i = 1:5
imagesc(fin_res(1).dat(i,2:end))
axis image off
pause
end



%% response
tmp = fin_res(3).resp(18).windows;
tmp(tmp > 10) = 0;
plot(fin_res(3).resp(18).x(100:end), tmp(:,100:end));
xlim([-20 30]);



%% response natural
%load '/Volumes/LAB/vega/EMGdat/EMGv010507.mat'
%saveas(gcf, ['~/Documents/uni/yifat_lab/results/poster/natural_response/' int2str(k) int2str(j) '.eps']);

files = what('/Volumes/LAB/vega/EMGdat/');
for k = 1:length(files.mat)
   load(['/Volumes/LAB/vega/EMGdat/' files.mat{k}]);
for j = 1:8
    pro = [];
    for i = 1:length(emgpsth)
        pro = [pro; emgpsth(i).hand(1).target(j,:)];
    end
    plot(pro');
    pause
end
end
  



