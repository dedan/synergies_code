function sep_results = sep_results_handposition(res_info, res_matrix)


pro_count = 1;
sup_count = 1;
mid_count = 1;
sep_results = struct;

for i = 1:length(res_info)
   switch(res_info(i).hand)
      case{1}
         sep_results(1).dat(pro_count,:)  = res_matrix(i,:);
         sep_results(1).info(pro_count)   = res_info(i).dat;
         pro_count = pro_count +1;
      case{2}
         sep_results(2).dat(sup_count,:)  = res_matrix(i,:);
         sep_results(2).info(sup_count)   = res_info(i).dat;
         sup_count = sup_count +1;
%       case{3}
%          sep_results(3).dat(mid_count,:)  = res_matrix(i,:);
%          sep_results(3).info(mid_count)   = res_info(i).dat;
%          mid_count = mid_count +1;
      otherwise
         disp(['now valid handpos, disregard: ' int2str(res_info(i).ses.id) ' -- ' int2str(res_info(i).ses.subsession) ' for further analysis']);
   end
end
