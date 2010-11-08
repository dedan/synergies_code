function sep_results = sep_results_handposition(resp)

% this function separated the responses according to their handposition
% (pronation and supination) I do not look for midsize here, as only a
% small number of sessions for vega has been recorded in midposition.
% NOTE for darma this has to be changed as we have only very little data

pro_count = 1;
sup_count = 1;
sep_results = struct;

for i = 1:length(resp)
   switch(resp(i).hand)
      case{1}
         sep_results(1).dat(pro_count,:)    = resp(i).response;
         sep_results(1).dat_y(pro_count,:)  = resp(i).response_y;
         sep_results(1).resp(pro_count)     = resp(i);
         pro_count = pro_count +1;
      case{2}
         sep_results(2).dat(sup_count,:)    = resp(i).response;
         sep_results(2).dat_y(sup_count,:)  = resp(i).response_y;
         sep_results(2).resp(sup_count)     = resp(i);
         sup_count = sup_count +1;
      otherwise
         disp(['now valid handpos, disregard: ' int2str(resp(i).info.ses.id) ' -- ' int2str(resp(i).info.ses.subsession) ' for further analysis']);
   end
end
