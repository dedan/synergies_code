

% don't remember what I used this for


pd_folder = '/Volumes/NEW_DISK/EMGdat_pd/';

dir_list = dir([pd_folder 'EMG_*_v*01ee_ch*_hand*.mat']);

vms = [];

for tmp_dir = dir_list'
   parsed = regexp(tmp_dir.name, '.*?_ch(.*?)_hand(.*?)\.mat', 'tokens');
   ch     = eval(char(parsed{1}(1)));
   hand   = eval(char(parsed{1}(2)));
   tmp = load(tmp_dir.name, 'vm');
   vms = [vms; [ch hand tmp.vm]];
end

for i = 1:16
   act_chan = vms(vms(:,1) == i,:);
   pro_dc = act_chan(act_chan(:,2) == 1,3);
   sup_dc = act_chan(act_chan(:,2) == 2,3);
   pro_g  = act_chan(act_chan(:,2) == 1,4);
   sup_g  = act_chan(act_chan(:,2) == 2,4);
   figure(1);
   subplot(4,4,i);
   plot(pro_dc);
   figure(2);
   subplot(4,4,i);
   plot(sup_dc);
   
   figure(3);
   subplot(4,4,i);
   plot(pro_g);
   figure(4);
   subplot(4,4,i);
   plot(sup_g);
end
