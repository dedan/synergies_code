function allall = nonstim_syns_svd()


% ich muss alle in eine matrix packen, die muskeln nach meiner schablone
% aussuchen, erst nach dem remapping (erst nach und ohne 250207)

Nitr = 100;
Ns = 3;
%NOTE config
e2take = find(logical([0  1  1  1  1  1  1  1  1  1  1  1  0  0  0  0 ]));


            % NOTE, bin mir nicht sicher ob ich nur die aussagekrŠfiten
            % matrizen will, wenn ich die svd dann spŠter fŸr alle mache

only_sig = 0;

res = struct('Wopt', [], 'Woptmax', struct, 'not_taken', []);

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
            R = S(1) / sum(diag(S));
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



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Rdat, Rmax, Wopt, Woptmax] = extract_syng( mat, mathand, Ns, Nitr)

Score1 = 99999;
Score2 = [99999 99999];
Wopt = [];
Woptmax = struct;

for z1 = 1:Nitr,
   [W,H,score1]=nmf(mat', Ns);
   if score1 < Score1
      Wopt  = W';
      Score1 = score1;
   end
   for z=1:length(mathand),
      currmat = mathand(z).data;
      H1 = extract_coeff( currmat', W,Ns);
      Rdat(z1,z) = get_resid( currmat, W, H1);

      % compraing the W taken from the 16 targets with W1 taking
      % using 8 targets only
      [Wmax,Hmax,score2]=nmf(currmat', Ns);
      Rmax(z1,z) = get_resid( currmat, Wmax, Hmax);
      if score2 < Score2(z),
         Woptmax(z).W = Wmax';
         Score2(z) = score2;
      end;
   end
end


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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function r = get_resid( matin, w, h)

SS = sum(sum(matin));
D = matin' - w*h;
r = sum(sum(abs(D)))/SS;
r= (1-r)*100;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function  Wout = average_over_sim( W1, W2)

L = size(W1,2);
for i=1:L,
   for j=1:L,
      cout(i,j) = pdist([(W1(:,i)) (W2(:,j))]');
   end
end

for i=1:L,
   [i,j]=find(cout == min(min(cout)));
   indx2(i)=j;
   cout(:,j)=99;
   cout(i,:)=99;
end

if length(unique(indx2))~= L,
   disp('Serious problem in matching');
end
Wout = (W1+W2(:,indx2))/2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Rmat, Rmathand] = get_random_matrices( mat, mathand)

Rmat = PermuteMatrix( mat );
for i=1:length(mathand),
   mat1 = mathand(i).data;
   Rmathand(i).data=PermuteMatrix(mat1);
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function omat = PermuteMatrix( mat );
omat = reshape(mat,1,[]);
indx = rand(size(omat));
[tmp,ri]=sort(indx);
omat = omat(ri);
omat=reshape(omat,size(mat));


