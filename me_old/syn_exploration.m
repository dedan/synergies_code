function [errors grouped_syns] = syn_exploration(group, dim, iter)

errors = zeros(1,dim);

% do the iterations
syns = struct;
for i = 1:iter
   syns(i).dat = NMF(group',dim)';
   syns(i).not_taken = 1:dim;
end

% for each synergist group (dimension) write the first entry, here I chose
% just to use one vector from the first iteration and group the others
% according to them. I am aware that there is a problem if these are
% outliers. Use the Neural gas version of this function !
grouped_syns = struct;
for i = 1:dim
   grouped_syns(i).dat(1,:) = syns(1).dat(i,:);
end


% such für jede dimension den nächstähnlichen, noch nicht vergebenen vector
% in allen anderen dimensionen
for i = 1:dim
   tmp = syns(1).dat(i,:);
   dif = 100000000000000;
   pos = -1;
   for j = 2:iter
      for k = find(syns(j).not_taken)
         cmp_dif = sum(tmp-syns(j).dat(k,:));
         if(cmp_dif < dif)
            dif = cmp_dif;
            pos = k;
         end
      end
      grouped_syns(i).dat(j,:) = syns(j).dat(pos,:);
      syns(j).not_taken(pos) = 0;
   end
end
      
for i = 1:dim
   errors(i) = error_function(mean(grouped_syns(i).dat), grouped_syns(i).dat);
end
