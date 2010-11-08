function showALLsyns_rnd()


% ich muss alle in eine matrix packen, die muskeln nach meiner schablone
% aussuchen, erst nach dem remapping (erst nach und ohne 250207)

%NOTE config
config = struct;

Nitr = 100;
Ns = 3;
e2take = find(logical([0  1  1  1  1  1  1  1  1  1  1  1  0  0  0  0 ]));
config.n_used_channels = length(e2take);

path = '/Volumes/lab/vega/EMGdat/';
outpath = '~/Documents/uni/yifat_lab/results/nonevoked_syns/';
% path = 'F:\vega\EMGdat\';
% outpath='F:\out\';

res = struct('runs', struct, 'Wopt', [], 'Woptmax1', [], 'Woptmax2', [], 'Wopt_rnd', [], 'data', struct);
fdat = dir([path 'EMGv*.mat']);
thresh = 70;


fdat  = sortfiles(fdat);
rc    = 1;

for i=1:length(fdat),

   
   % NOTE only for debug
   if i < 25
      continue
   end
   
   data = load([path char(fdat(i))]);
   
   % only when data for different handpositions available
   if length(data.chdata) > 1

      % get the data in mat that contains all and in a struct mathand
      % separated by handposition
      [mat, mathand] = get_matrices( data.chdata, e2take);
      if ~isempty(mat)
         
         % random matrices for comparison
         [Rmat, Rmathand]  = get_random_matrices( mat, mathand);

         % put in struct what we got so far
         res(rc).data.mat  = mat;
         res(rc).data.name = char(fdat(i));
         res(rc).data.sig  = -1;

         % explained variance tests
         res(rc).data.r_nmf         = test_resid_nmf( mat,Nitr );
         res(rc).data.r_rnd_nmf     = test_resid_nmf( Rmat,Nitr );
         [U S]                      = svd(mat);
         res(rc).data.eigen         = diag(S);
         [U S]                      = svd(mat);
         res(rc).data.eigen_rnd     = diag(S);
         res(rc).data.r_pca         = test_resid_pcaica(mat, config);
         res(rc).data.r_rnd_pca     = test_resid_pcaica(Rmat, config);
         res(rc).data.r_svd         = test_resid_svd(mat, config);
         res(rc).data.r_rnd_svd     = test_resid_svd(Rmat, config);
         
         
         
            
            subplot 221
            ax = 1:length(R);
            cum = zeros(length(res(rc).data.eigen),1);
            for k = 1:length(res(rc).data.eigen)
               cum(k) = sum(res(rc).data.eigen(1:k));
            end
            Rsvd = [( 100 - (cum / sum(diag(S))) * 100)' zeros(1,length(R)-length(diag(S)))];
            plot(ax,R,'.b',ax,Rr,'.r',ax,Rsvd,'.g');
            title('Variance-dimension dependency');
            xlabel('Dimensions');
            ylabel('Explained variance');

            % compute it five times to see stability of nmf results
            for j = 1:5
               [Rdat, Rmax, Wopt, Woptmax] = extract_syng( mat, mathand, Ns, Nitr);
               [Rrdat, Rrmax, Wopt_rnd, Woptmax_rnd]  = extract_syng(Rmat, Rmathand, Ns, Nitr);
               if j == 1
                  res(rc).Wopt     = Wopt;
                  res(rc).Woptmax1 = Woptmax(1).W;
                  res(rc).Woptmax2 = Woptmax(2).W;
                  res(rc).Wopt_rnd = Wopt_rnd;
               end
               res(rc).runs(j).Wopt      = Wopt;
               res(rc).runs(j).Woptmax1  = Woptmax(1).W;
               res(rc).runs(j).Woptmax2  = Woptmax(2).W;               
               res(rc).runs(j).Wopt_rnd      = Wopt_rnd;
               res(rc).runs(j).Woptmax1_rnd  = Woptmax_rnd(1).W;
               res(rc).runs(j).Woptmax2_rnd  = Woptmax_rnd(2).W;               
            end

            subplot 222
            %NOTE schauen warum ich bei Rmax nicht zwei balken bekomme
            val = [mean(Rdat);mean(Rmax); mean(Rrdat); mean(Rrmax)]';
            eval = [std(Rdat);std(Rmax); std(Rrdat); std(Rrmax)]';
            barweb( val,eval,.5,[], [],'Hand Position','Explained variance', 'default', [], {'dat','max','Rdat', 'Rmax'})

            subplot 223
            bar(Wopt'),
            xlabel('Muscle');
            ylabel('Activation');
            legend('Syn-1','Syn-2','Syn-3');

            subplot 224
            plot(mat,'-+');
            title('Muscle tuning');
            xlabel('Targets');
            ylabel('Muscle Activation');

            figname = char(fdat(i));
            pos = min(findstr(figname,'.'));
            %NOTE config
            figname = [outpath 'OUT' figname(4:pos) 'fig.tif'];
            saveas(gcf, figname, 'tiff');

            clear  Rmax Rdat Rrdat Rrmax
   
         end
         rc = rc +1;
      end
   end
end

save([outpath 'all_syns'], 'res');



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
   %NOTE no midpostion
   for z=1:2
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
for k=1:length(chdata),
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


