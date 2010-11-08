function [agg agall] = get_nonstim_data()


% ich muss alle in eine matrix packen, die muskeln nach meiner schablone
% aussuchen, erst nach dem remapping (erst nach und ohne 250207)

%NOTE config
e2take = find(logical([0  1  1  1  1  1  1  1  1  1  1  1  0  0  0  0 ]));


% NOTE, bin mir nicht sicher ob ich nur die aussagekrŠfiten
% matrizen will, wenn ich die svd dann spŠter fŸr alle mache

only_sig = 0;

fdat = dir(['/Volumes/lab/vega/EMGdat' filesep 'EMGv*.mat']);
thresh = 70;


fdat = sortfiles(fdat);
res_count = 1;
agg = struct;
agall = struct;

for i=1:length(fdat),

   data = load(['/Volumes/lab/vega/EMGdat' filesep char(fdat(i))]);

   % only when data for different handpositions available
   if length(data.chdata) > 1

      % get the data in mat that contains all and in a struct mathand
      % separated by handposition
      [mat, mathand] = get_matrices( data.chdata, e2take);

      if ~isempty(mat)

         if only_sig

            [U S V] = svd(mat);
            R = (S(1) / sum(diag(S))) * 100;
            disp(['Rank1 explains --> ' num2str(R(1))]);
            if R(1) < thresh,  % threshold is 80 or darma and 75 for Vega
               disp('Rank-1 of the EMG matrix is insufficient to explain the data');

               agg(res_count).both    = mat;
               agg(res_count).pro     = mathand(1).W;
               agg(res_count).sup      = mathand(2).W;
               agall.pro = [agall.pro; mathand(1).W];
               agall.sup = [agall.sup; mathand(2).W];
               agall.both = [agall.both; mat];
               res_count = res_count +1;

            end
         else
            agg(i).both    = mat;
            agg(i).pro     = mathand(1).W;
            agg(i).sup      = mathand(2).W;
            agall.pro = [agall.pro; mathand(1).W];
            agall.sup = [agall.sup; mathand(2).W];
            agall.both = [agall.both; mat];
         end
      end
   end
end

%save('~/Documents/uni/yifat_lab/results/nonevoked_syns/all_syns', 'res');





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function odir = sortfiles( idir )

for i=1:length(idir),
   curname = char(idir(i).name);
   odir(i) = {curname};
   curname = curname(4:end);
   Nday = str2num( curname(2:3));
   Nmonth = str2num(curname(4:5));
   Ndates(i,1) = Nmonth;
   Ndates(i,2) = Nday;
end
[indx,ilist]=sortrows(Ndates);
odir = (odir(ilist));



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [mat, mathand] = get_matrices( chdata, e2take)

mat =[];
mathand = [];

if length(e2take) < 5 ,% we are using at least 8 muscles
   disp('Not enought EMGs for analysis');
   return;
end

% only pronation and supination
for k=1:2
   mat2add = chdata(k).mat';
   if size(chdata(k).mat,1) > length(e2take)
      mat2add = mat2add(:,e2take);
   end
   mathand(k).data = mat2add;
   mat = [mat' mat2add']';
end

