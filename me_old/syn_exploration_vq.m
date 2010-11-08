function expl_res = syn_exploration_vq(group, config)

dim            = config.dimensions;
iter           = config.Niter_exploration;

expl_res                = struct;
expl_res.all_patterns   = zeros(dim*iter,config.n_used_channels);

% führe die iterationen durch
syns = struct;
for i = 1:iter
   syns(i).dat = NMF(group',dim)';
   syns(i).not_taken = 1:dim;
   expl_res.all_patterns(dim*i - dim +1:dim*i,:) = syns(i).dat;
end


% find prototypes by k_means clustering
[IDX,protos]      = kmeans(expl_res.all_patterns,dim);
expl_res.protos   = protos;
expl_res.idx      = IDX;

grouped_syns = struct;

for i = 1:dim
   tmp = protos(i,:);
   for j = 1:iter
      dif = 1000000000000000000000;
      pos = -1;
      for k = find(syns(j).not_taken)
         cmp_dif = pdist([tmp; syns(j).dat(k,:)]);
         if(cmp_dif < dif)
            dif = cmp_dif;
            pos = k;
         end
      end
      grouped_syns(i).dat(j,:) = syns(j).dat(pos,:);
      syns(j).not_taken(pos) = 0;
   end
end

expl_res.grouped_syns = grouped_syns;
