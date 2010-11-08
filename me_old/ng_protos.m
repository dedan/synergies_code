function protos = ng_protos(all_patterns, dim, eps, lam, iter, eps_d, lam_d)


protos       = rand(dim,size(all_patterns,2)) * max(max(all_patterns));
dists        = zeros(1,dim);
for k = 1:iter
   for i = 1:size(all_patterns,1)
      for j = 1:dim
         dists(j) = norm(protos(j,:) - all_patterns(i,:))^2;
      end
      [sorted_d idx] = sort(dists);
      for j = 1:dim
         protos(j,:) = protos(j,:) + eps * exp( -(idx(j)-1) / lam) * (all_patterns(i,:) - protos(j,:));
      end
   end
   eps = eps * eps_d;
   lam = lam * lam_d;
end