function grouped = group(data, fieldname)

debug_info  = 0;
n_iter      = 100;
big_number  = 10000000000000;
group_size  = size(data(1).(fieldname),1);
grouped     = struct;
add_info    = struct;
flat        = [];

for i = 1:length(data)
   flat                    = [flat; data(i).(fieldname)]; %#ok<AGROW>
   add_info(i).not_taken   = true(1,group_size);
end

dif = big_number;
for i = 1:n_iter
   % NOTE verschiedene emptyactions testen
   [tmp_idx tmp_c err] = kmeans(flat, group_size,'emptyaction', 'singleton', 'display', 'off');
   if debug_info
      bla = histc(tmp_idx,1:3);
      disp(['err: ' num2str(sum(err)) ' dist: ' int2str(bla(1)) '_' int2str(bla(2)) '_' int2str(bla(3))]);
   end
   if sum(err) < dif
      dif   = sum(err);
      idx   = tmp_idx;
      c     = tmp_c;
   end
end

grouped.idx    = idx;
grouped.center = c;
if debug_info
   bla = histc(idx,1:3);
   disp(['final err: ' num2str(sum(dif)) ' dist: ' int2str(bla(1)) '_' int2str(bla(2)) '_' int2str(bla(3))]);
end
for i = 1:group_size
   for j = 1:length(data)
      dif = big_number;
      pos = -1;
      for k = find(add_info(j).not_taken)
         cmp_dif = pdist([c(i,:); data(j).(fieldname)(k,:)]);
         if(cmp_dif < dif)
            dif = cmp_dif;
            pos = k;
         end
      end
      grouped(i).dat(j,:)        = data(j).(fieldname)(pos,:);
      add_info(j).not_taken(pos) = false;
   end
end
