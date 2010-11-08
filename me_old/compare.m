
% compare all the results

path        = '~/Documents/uni/yifat_lab/results/data/';
files       = what(path);
n_files     = length(files.mat);
n_first     = 3;
figure(1);

for k = 1:2
   figure(k);
   for i = 1:n_files;
      dat = load([path files.mat{i}]);
      if k == 1
         disp(dat.fin_res(1).config.comment);
      end
      for j = 1:n_first
         subplot(n_files, n_first+1, ((i-1) * (n_first+1)) +j);
         bar(dat.fin_res(k).protos(dat.fin_res(k).syn_order(j),:));
      end
      subplot(n_files, n_first+1, ((i-1) * (n_first+1)) +4);
      imagesc(dat.fin_res(k).protos);
   end
end