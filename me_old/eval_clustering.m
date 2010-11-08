

% this script is to evaluate the stability of the prototyping algorith.
% runs the clustering several times with the same amount of iterations.
% afterwards compares the computed synergists

addpath(genpath('..')); %#ok<MCAP>

hand_pos          = 1;
dimensions        = 6;
Niter_exploration = 1250;
Niter_repeat      = 10;
N_best            = 3;

fin_res           = struct;
best_syns         = zeros(N_best * Niter_repeat, 16);
agg               = [];

j = hand_pos;
% do the clustering several times
for k = 1:Niter_repeat

   % compute the synergists for the final results and also group it
   % randomly to get compare the sorting with
   fin_res(j).dat   = sep_results(j).dat(remap_idx(j).after,:);
   fin_res(j).dim   = sep_results(j).dim_after;
   [protos grouped_syns_vq]   = syn_exploration_vq(fin_res(j).dat, dimensions, Niter_exploration, 0.5, 10, 500, 0.995, 0.99);
   fin_res(j).protos          = protos;
   fin_res(j).grp_syns        = grouped_syns_vq;

   % evaluate the errors
   for i = 1:dimensions
      fin_res(j).mean_syns(i,:)       = mean(grouped_syns_vq(i).dat);
      fin_res(j).errors(i)            = error_function(protos(i,:), grouped_syns_vq(i).dat);
   end

   [sorted_errors_vq  idx_vq]    = sort(fin_res(hand_pos).errors);
   sorted_syns                   = fin_res(j).mean_syns(idx_vq,:);
   sorted_protos                 = fin_res(j).protos(idx_vq,:);

   agg   = [agg; sorted_syns(1:N_best,:)];
   
end

figure(1);
imagesc(agg);
