
% analysis of the synergist found in a non-stimulation condition

%% analysis of residual tests

if 1

load ~/Documents/uni/yifat_lab/results/nonevoked_syns/all_syns.mat

res_fil = struct('runs', struct, 'Wopt', [], 'Wopt_rnd', [], 'data', struct,'Woptmax1', [], 'Woptmax2', []);

nmf_svd_all = [];
nmf_svd_sig = [];
nmf_svd_non = [];
rank1       = [];
all_eigen   = [];
all_eigen_sig = [];

j = 1;
for i = 1:length(res)
   all_eigen = [all_eigen res(i).data.eigen / sum(res(i).data.eigen)];
   r = res(i).data.resid(1:11);
   cum = zeros(length(res(i).data.eigen),1);
   for k = 1:length(res(i).data.eigen)
      cum(k) = sum(res(i).data.eigen(1:k));
   end
   e = (100 - (cum / sum(res(i).data.eigen)) * 100)';
   nmf_svd_all = [nmf_svd_all [r; e]];

   rank1 = [rank1 [r(1); (100 - (cum(1) / cum(end)) * 100)' ]];

   if res(i).data.sig
      all_eigen_sig = [all_eigen_sig res(i).data.eigen / sum(res(i).data.eigen)];
      nmf_svd_sig = [nmf_svd_sig [r; e]];
      res_fil(j)  = res(i);
      j = j+1;
   else
      nmf_svd_non = [nmf_svd_non [r; e]];
   end
end

figure
subplot 221
plot(nmf_svd_all(1,:), nmf_svd_all(2,:),'.');
axis square
axis([0 100 0 100])
xlabel('nmf');
ylabel('svd');
title('remaining errors - all')

subplot 222
plot(rank1(1,:), rank1(2,:),'.');
axis square
axis([0 100 0 100])
xlabel('nmf');
ylabel('svd');
title('remaining errors - rank1')

subplot 223
plot(nmf_svd_sig(1,:), nmf_svd_sig(2,:),'.');
axis square
axis([0 100 0 100])
xlabel('nmf');
ylabel('svd');
title('remaining errors - only rank1 < thres')

subplot 224
plot(nmf_svd_non(1,:), nmf_svd_non(2,:),'.');
axis square
axis([0 100 0 100])
xlabel('nmf');
ylabel('svd');
title('remaining errors - only rank1 > thres')

figure
subplot 211;
plot(all_eigen);
title('eigenvalues');

subplot 212;
plot(all_eigen_sig);
title('eigenvalues - only rank1 < thres');


res = res_fil;


end

%% stability of nmf graphically

% all_nmf     = [];
% all_nmf_rnd = [];
% centers     = struct;
%
% % soll zeigen, dass nmf stabiler, wenn struktur in daten und nicht random
%
% dist_c = 1;
% for i = 1:length(res)
%    grouped     = group(res(i).runs, 'Wopt');
%    grouped_rnd = group(res(i).runs, 'Wopt_rnd');
%    centers(i).syn   = grouped(1).center;
%    centers(i).rnd   = grouped_rnd(1).center;
%    for j = 1:length(grouped)
%       all_nmf        = [all_nmf; grouped(j).dat ];
%       group_avg_dist(dist_c) = mean(pdist(grouped(j).dat));
%
%       all_nmf_rnd        = [all_nmf_rnd; grouped_rnd(j).dat ];
%       group_avg_dist_rnd(dist_c) = mean(pdist(grouped_rnd(j).dat));
%       dist_c = dist_c +1;
%    end
% end
%
% figure
% imagesc(all_nmf);
% figure
% imagesc(all_nmf_rnd);
% figure;
% subplot 511
% imagesc(all_nmf);
% subplot 512
% hist(group_avg_dist,10);
% subplot 513
% imagesc(all_nmf_rnd);
% subplot 514
% hist(group_avg_dist_rnd,10);
%
%
% % jetzt avg_dist zwischen den prototypen für verschiedene sessions
%
% grouped = group(centers, 'syn');
% grouped_rnd = group(centers, 'rnd');
%
% for i = 1:length(grouped)
%    syn_avg_dist(i) = mean(pdist(grouped(i).dat));
%    rnd_avg_dist(i) = mean(pdist(grouped_rnd(i).dat));
% end
%
% subplot 515
% bar([syn_avg_dist; rnd_avg_dist]);
%
%
%
% %% consistency over sessions
%
% % da es so aussieht als sei die nmf stabil, nehme ich jetzt einfach das
% % erste der 5 ergebnisse
%
% grouped = group(res, 'Wopt');
% flat = [];
% for i = 1:length(grouped)
%    flat = [flat; grouped(i).dat];
% end
%
% figure;
% imagesc(flat);


%% ich mache jetzt erst mal mit den protypen weiter

% ich muss jetzt auf den unterschied der projection auf den raum von c
% und den orthogonal complementären raum von c schauen.

% c        = grouped.center;
% basis    = orth(c');
% ortcmp   = null(basis');
% p        = basis * basis';
% p2       = ortcmp * ortcmp';
%
% load ~/Documents/uni/yifat_lab/results/30122008__213/data.mat
%
% figure;
% for i = 1:length(fin_res)
%    dat = fin_res(i).dat';
%
%    sep_grouped = group(res, ['Woptmax' int2str(i)]);
%    sep_c       = sep_grouped.center;
%    sep_basis   = orth(sep_c');
%    sep_p       = sep_basis * sep_basis';
%    sep_ortcmp  = null(sep_basis');
%    sep_p2      = sep_ortcmp * sep_ortcmp';
%
%
%    dat = reshape(dat(randperm(length(dat(:)))), size(dat) );
%    for j = 1:length(dat)
%       dat(:,j) = dat(:,j) ./ norm(dat(:,j));
%    end
%
%
%    on_b     = p * dat;
%    on_o     = p2 * dat;
%    on_sep_b = sep_p * dat;
%    on_sep_o = sep_p2 * dat;
%
% %    for j = 1:length(dat)
% %       dist_b(j)      = norm(on_b(:,j));
% %       dist_o(j)      = norm(on_o(:,j));
% %       dist_sep_b(j)  = norm(on_sep_b(:,j));
% %       dist_sep_o(j)  = norm(on_sep_o(:,j));
% %    end
% %
%    for j = 1:length(dat)
%       dist_b(j)      = sum(on_b(:,j));
%       dist_o(j)      = sum(on_o(:,j));
%       dist_sep_b(j)  = sum(on_sep_b(:,j));
%       dist_sep_o(j)  = sum(on_sep_o(:,j));
%    end
%
%
%
%    % NOTE hier auch zahlen, werte, nicht nur verteilungen ausgeben
%
%    subplot(4,2,i);
%    plot(dist_b,dist_o,'.');
%    axis square
%    axis([0 1 0 1]);
%    xlabel('distance on synergis basis');
%    ylabel('distance on orthogonal complement');
%
%    subplot(4,2,i+2);
%    plot(dist_sep_b,dist_sep_o,'.');
%    axis square
%    axis([0 1 0 1]);
%    xlabel('distance on synergis basis');
%    ylabel('distance on orthogonal complement');
%
%
%    subplot(4,2,4+i)
%    hist(dist_b,10);
%    subplot(4,2,6+i);
%    hist(dist_o,10);
% end
%


%% relations with bootstrapping

% load the evoked response data
load ~/Documents/uni/yifat_lab/results/30122008__213/data.mat
dat = fin_res(1).dat';
%dat = fin_res(1).protos';

dat = fin_res(1).pca_syns';

n_iter = 1000;
dist_length = 100000;


% create basis and projections
grouped  = group(res, 'Wopt');
c        = grouped.center;
basis    = orth(c');
ortcmp   = null(basis');
space    = [basis ortcmp];
p        = space \ (basis * basis');
p2       = space \ (ortcmp * ortcmp');

data_ratio  = sum(sum(abs(p * dat))) / sum(sum(abs(p2 * dat)));
shuf_ratios = NaN(1,n_iter);
rand_ratios = NaN(1,n_iter);
test_ratios = NaN(1,n_iter);

% bootstrapping loop
for i = 1:n_iter

   %create test_space and compute ratio
   x        = linsolve(basis,dat);
   mean_x   = mean(x,2);
   std_x    = std(x,0,2);
   x_basis  = NaN(size(x));
   for j = 1:size(x,1)
      distribution = (randn(1,dist_length) * std_x(j)) + mean_x(j);
      perms        = randperm(dist_length);
      x_basis(j,:) = distribution(perms(1:size(x,2)));
   end
   dat_test       = basis * x_basis;
   % put noise on test data, otherwise very high value
   noise_10p_gaus = randn(size(dat_test)) * mean(dat_test(:)) * 0.5;
   dat_test       = dat_test + noise_10p_gaus;
   test_ratios(i) = sum(sum(abs(p * dat_test))) / sum(sum(abs(p2 * dat_test)));

   %shuffle data and compute ratio
   dat_shuf       = reshape(dat(randperm(length(dat(:)))), size(dat) );
   shuf_ratios(i) = sum(sum(abs(p * dat_shuf))) / sum(sum(abs(p2 * dat_shuf)));
   
   %create random data and compute ratio
   dat_rand       = rand(size(dat));
   rand_ratios(i) = sum(sum(abs(p * dat_rand))) / sum(sum(abs(p2 * dat_rand)));
end


% estimate offset for my basis
n     = 1000;
dim   = fin_res(1).config.n_used_channels;
rel_res = zeros(3,floor(dim/2));
for sub_dim = 1:floor(dim/2)
   space    = [orth(c') null(basis')];
   part  = space(:,1:sub_dim);
   comp  = null(part');
   p = space \ (part  * part');
   p2 = space \ (comp  * comp');
   pos = rand(dim,n);
   rel_res(1,sub_dim) = sum(sum(abs(p * pos))) / sum(sum(abs(p2 * pos)));
   rel_res(2,sub_dim) = sub_dim / (dim - sub_dim);
   rel_res(3,sub_dim) = rel_res(1,sub_dim) - rel_res(2,sub_dim);
end

disp(['b zu o normal: ' num2str(data_ratio)]);
disp(['mean ratio from shuffle: ' num2str(mean(shuf_ratios))]);
disp(['mean rand ratio: ' num2str(mean(rand_ratios))]);
disp(['mean test ratio: ' num2str(mean(test_ratios))]);

disp(' ');
disp(['dim(b) zu dim(o): ' num2str(rank(basis) / rank(ortcmp))]);
disp(['mean rand ratio minus estimated offset: ' num2str(mean(rand_ratios) - rel_res(3,size(basis,2)))]);
disp(['mean shuf ratio minus estimated offset: ' num2str(mean(shuf_ratios) - rel_res(3,size(basis,2)))]);
disp(['data ratio minus estimated offset: ' num2str(data_ratio - rel_res(3,size(basis,2)))]);
disp(['mean test ratio minus estimated offset: ' num2str(mean(test_ratios) - rel_res(3,size(basis,2)))]);



figure(8);
subplot 511
hist(shuf_ratios);
subplot 512
hist(rand_ratios);
subplot 513
hist(test_ratios);
subplot(5,1,4:5);
plot(rel_res');
hold on
plot(ones(length(rel_res),1) * rel_res(3,size(basis,2)),'--');
hold off









% %% synergie pd plot
% % plot the synergies over all sessions for both handpositions
% figure;
%
% for j = 1:length(fin_res)
%    sep_grouped = group(res, ['Woptmax' int2str(j)]);
%    c           = sep_grouped.center;
%    n_syns = size(c,1);
%    for k = 1:n_syns
%       subplot(2,n_syns,k + (j-1)*3);
%       rose_agg = [];
%       used_pds = fin_res(j).channel_pd(fin_res(j).config.channels2take);
%       for i = 1:fin_res(j).config.n_used_channels
%          rose_agg = [rose_agg ones(1,floor(c(k,i) * 100)) * used_pds(i)]; %#ok<AGROW>
%       end
%       h_fake = rose(ones(1,250));
%       hold on;
%       rose(rose_agg);
%       set(h_fake, 'Visible', 'Off');
%       title(['pos: ' int2str(j) ' syn: ' int2str(k)]);
%    end
% end
%
%
% %% check pds on consistency between handpositions
%
% % figure
% % for j = 1:length(fin_res)
% %    subplot(length(fin_res),1,j);
% %    rose_agg = [];
% %    used_pds = fin_res(j).channel_pd(fin_res(j).config.channels2take);
% %    for i = 1:fin_res(j).config.n_used_channels
% %       rose_agg = [rose_agg 100 * used_pds(i)]; %#ok<AGROW>
% %    end
% %    rose(rose_agg);
% % end
%
%
% figure
% for j = 1:length(fin_res)
%    subplot(length(fin_res),1,j);
%    used_pds = fin_res(j).channel_pd(fin_res(j).config.channels2take);
%    [x y] =  pol2cart(used_pds, ones(size(used_pds)));
%    compass(x, y);
% end
%

