function [W,H, score] = NMF(V,k)
% INPUT: V is m*n non-negative  matrix , k is the reduced dimension 
% OUTPUT: W is m*k,  H is k*n non-negative matrices  s.t. V ~ WH
[m,n]  = size(V);
W = rand(m,k);
H = rand(k,n);

maxiter =50 ;

for i= 1:maxiter
    H = H.*(W'*V)./(W' *W * H  + 1e-10)  ;
    W = W.* (V * H')  ./ ( W * H * H' + 1e-10) ; 
end


dif = V - W*H;
score = sum(sum(sqrt(dif.*dif)));
