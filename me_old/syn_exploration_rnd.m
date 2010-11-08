function [grouped_syns] = syn_exploration_rnd(group, config)

dim   = config.dimensions;
iter  = config.Niter_exploration;

% führe die iterationen durch
syns = struct;
for i = 1:iter
   syns(i).dat = NMF(group',dim)';
   syns(i).not_taken = 1:dim;
end


% such für jede dimension den nächstähnlichen, noch nicht vergebenen vector
% in allen anderen dimensionen
for i = 1:dim
   for j = 1:iter
      not_takens = find(syns(j).not_taken);
      rand_ints = randperm(length(not_takens));
      pos = not_takens(rand_ints(1));
      grouped_syns(i).dat(j,:) = syns(j).dat(pos,:);
      syns(j).not_taken(pos) = 0;
   end
end
