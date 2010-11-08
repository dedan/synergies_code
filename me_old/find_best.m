

%% back nmf

V = all';
W = prot';
k = size(prot,1);

[m,n]  = size(V);
H = rand(k,n);

maxiter =50 ;

for i= 1:maxiter
    H = H.*(W'*V)./(W' *W * H  + 1e-10)  ;
    W = W.* (V * H')  ./ ( W * H * H' + 1e-10) ; 
end


dif = V - W*H;
score = sum(sum(sqrt(dif.*dif)))

dif_w = W - prot';
score_w = sum(sum(sqrt(dif_w.*dif_w)))




%% analysis

res = ones(k,1);

for i = 1:k
   choose = logical(find(1:11));
   choose(i) = 0;
   dif = sum(sum(abs(V - W(:,choose) * H(choose,:))));
   disp(['prototype: ' int2str(i) ' fehler: ' num2str(dif)]);
   res(i) = dif;
end

[sorted idx] = sort(res, 'descend');

figure(5);
plot(sorted);

err = zeros(k,1);

choose = logical(find(1:k));
for i = 1:k
   choose(idx(i)) = 0;
   choose = logical(choose);
   err(i) = sum(sum(abs(V - W(:,choose) * H(choose,:))));
end

figure(6);
plot(err);





