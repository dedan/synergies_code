
%% check how stable the idea of usin NMF backwards is

res = [];
[w h score] = NMF(fin_res(1).dat', 10);

for i = 1:100
   [w1 h1 score1] = NMF(fin_res(1).dat', 10, w);
   res = [res sum(sum(abs(h-h1)))];
   score = [score score1]; %#ok<AGROW>
end

figure(1);
subplot 121;
plot(score);
subplot 122;
plot(res);



%% check the stability of k-means
% and how the error and the distribution of cluster size correlates

res_err = [];
res_var = [];
old = 1000000000000000;
min = [];

for i = 1:100
   [IDX C sumd] = kmeans(expl_res(1).all_patterns, 10);
   res_err = [res_err sum(sumd)];
   if sum(sumd) < old
      old = sum(sumd);
      min = histc(IDX,1:10);
   end
   res_var = [res_var std(histc(IDX,1:10))];
end

figure(2)
subplot 121;
plot(res_err);
subplot 122;
plot(res_err,res_var,'.');
min


%% playing with linkage

y = pdist(expl_res(1).all_patterns);
z = linkage(y,'complete');
dendrogram(z,'colorthreshold','default');


%% split kmeans
 options = [1 0.001 0.001 0 1]; 
 options(14)=1000;
 
[centers, options, post,errlog] = kmeans_split(expl_res(1).all_patterns, options, 10, 'dist2');
size(centers)
