
% look at the consistency of the preferred directions over sessions


load ../../results/pd/allpd.mat

channels=[11 12 13 14 21 22 23 24 31 32 33 34];

figure(1);
title('pronation: blue -- supination: red');
for i = 1:length(channels)
   subplot(6,2,i);
   pds_act_channel = alldout(find(alldout(:,1) == channels(i)),[2 4] );
   pro = pds_act_channel(pds_act_channel(:,2) == 1);
   sup = pds_act_channel(pds_act_channel(:,2) == 2);
   if(~isempty(pds_act_channel))
      [pro_x pro_y] = pol2cart(pro, ones(size(pro)));
      [sup_x sup_y] = pol2cart(sup, ones(size(sup)));
      hold on;
      feather(pro_x,pro_y, 'b');
      feather(sup_x,sup_y, 'r');
      hold off;
      title(['channel: ' int2str(i-1)]);
   end
end
