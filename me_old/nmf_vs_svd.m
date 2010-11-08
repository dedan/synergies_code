
% comparison of residuals when using svd or nmf.
% the data which is used hier (fin_res) are computed by
% the first cells of the synergists/syn_analysis.m script

% number of iterations to overcome local maxima of nmf
Nitr = 50;

matin = fin_res(1).dat;

SS = sum(sum(matin));

for i=1:20,
   
   % svd
   [U S V]     = svd(matin);
   S(i+1:end,i+1:end)  = 0;
   dif_svd     = matin - (U*S*V');
   err_svd(i)  = sum(sum(abs(dif_svd))); 
   
   % repeat it for nmf
    for j=1:Nitr,
        [w,h]     = NMF(matin',i);
        dif_nmf   = matin' - w*h;
        err_nmf(j,i)    = sum(sum(abs(dif_nmf)));
    end
end

%v = median(v)/SS;

err_nmf = min(err_nmf)/SS;
err_svd = err_svd / SS;

err_nmf = (ones(size(err_nmf))-err_nmf)*100;
err_svd = (ones(size(err_svd))-err_svd)*100;

figure(1);
hold on;
plot(err_nmf, 'r');
plot(err_svd, 'g');


% % v =mean(v) /SS;
% % v= (ones(size(v))-v)*100;
% % 

