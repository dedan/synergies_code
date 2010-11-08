function pds = get_pds(config)

% compute the average preferred directions over all sessions, I don't use
% this at the moment although the feather plot showed nice consistency over
% sessions, see: './pd_consistency'
% Maybe could be usefull again when looking for an orientation map on the
% cortex

pds = NaN(2, config.n_used_channels);

% read in all the PDs, separated by hand position
pd_files = what(config.pd_folder);
count    = length(pd_files.mat);
all_pd1  = NaN(config.n_used_channels, count);
all_pd2  = NaN(config.n_used_channels, count);
for i = 1:count
   pd_file = char(pd_files.mat(i));
   load([config.pd_folder pd_file]);
   for j = 1:length(VMout);
      if ~isempty(VMout(j).hand)
         if ~isempty(VMout(j).hand(1).PD) && VMout(j).hand(1).PD ~= 0
            all_pd1(j,i) = VMout(j).hand(1).PD;
         end
         if ~isempty(VMout(j).hand(2).PD) && VMout(j).hand(2).PD ~= 0
            all_pd2(j,i) = VMout(j).hand(2).PD;
         end
      end
   end
end

disp([int2str(sum(sum(isnan(all_pd1)))) ' pd values are missing for pronation']);
disp([int2str(sum(sum(isnan(all_pd2)))) ' pd values are missing for supination']);

% compute the channel PDs (remove non-valid PD and build vector sum) this
% is done in cartesian coordinates and reconverted to polar coordinates
for i = 1:config.n_used_channels
   pro            = all_pd1(i, 1 == isfinite(all_pd1(i,:)));
   sup            = all_pd2(i, 1 == isfinite(all_pd2(i,:)));
   [pro_x pro_y]  = pol2cart(pro, ones(size(pro)));
   [sup_x sup_y]  = pol2cart(sup, ones(size(sup)));
   pds(1,i)       = cart2pol(sum(pro_x),sum(pro_y));
   pds(2,i)       = cart2pol(sum(sup_x),sum(sup_y));
end
