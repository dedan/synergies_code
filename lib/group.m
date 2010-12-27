function grouped = group(data_orig, fieldname)

% data: a struct where the field with the name fieldname contains matrices
% that for example represent results of different runs of an algorithm. 

% I used it for checking the stability of the outcome of nnmf. As the
% resulting vectors are not ordered like for example in pca, I had to order
% and group them

n_iter      = 100;
big_number  = 10000000000000;
group_size  = size(data_orig(1).(fieldname),1);
grouped     = struct;
add_info    = struct;
flat        = [];

for i = 1:length(data_orig)
   data(i).dat             = normr(data_orig(i).(fieldname));     %#ok<AGROW>
   flat                    = [flat; data(i).dat]; %#ok<AGROW>
   add_info(i).not_taken   = true(1,group_size);
end


dif = big_number;

% run kmeans several times as it is known to be not the most stable
% algorithm
warning off stats:kmeans:EmptyCluster
for i = 1:n_iter
   [tmp_idx tmp_c err] = kmeans(flat, group_size,'emptyaction', 'singleton', 'display', 'off');   
   if sum(err) < dif
      dif   = sum(err);
      idx   = tmp_idx;
      c     = tmp_c;
   end
end
warning on stats:kmeans:EmptyCluster


grouped.idx    = idx;
grouped.center = c;

gs = 1:group_size;

% for all data eintries (runs)
for j = 1:length(data)
    
    % choose group randomly in order to prevent bias in groups
   for i = gs(randperm(length(gs)));
      dif = big_number;
      pos = -1;
      
      % look for the closest to prototype entry in this group which is
      % not already taken
      for k = find(add_info(j).not_taken)
         cmp_dif = pdist([c(i,:); data(j).dat(k,:)]);
         if(cmp_dif < dif)
            dif = cmp_dif;
            pos = k;
         end
      end
      
      grouped(i).dat(j,:)        = data_orig(j).(fieldname)(pos,:);
      add_info(j).not_taken(pos) = false;
   end
end



