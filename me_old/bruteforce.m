
% learn params

old_value = 1000000000000000;
disp('test');
for eps = 0.1:0.1:10
   for lam = 0.1:0.5:50
      for d_eps = 0.5:0.01:0.99
         for d_lam = 0.5:0.01:0.99
            [errors3 grouped_syns3] = syn_exploration_vq(sep_results(1).dat, sep_results(1).dim, 100, eps, lam, 200, d_eps, d_lam);
            if(sum(errors3) < old_value)
               old_value = sum(errors3);
               disp(['eps: ' num2str(eps) ' lam: ' num2str(lam) ' d_eps: ' num2str(d_eps) ' d_lam: ' num2str(d_lam) ' err: ' num2str(old_value)]);
            end
         end
      end
   end
end