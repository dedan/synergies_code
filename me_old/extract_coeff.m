function H = extract_coeff(V,W,k)
% INPUT: V is m*n non-negative  matrix , W is an m*k non negat matrix. 
% OUTPUT: H is k*n non-negative matrices  s.t. V ~ WH
% % % [m,n]  = size(V);
% % % [tmp,k]=size(W);
% % % H = rand(k,n);
% % % 
% % % H = (W'*V);
% % % 


[m,n]  = size(V);
% W = rand(m,k);
H = rand(k,n);

maxiter =50 ;

for i= 1:maxiter
    H = H.*(W'*V)./(W' *W * H  + 1e-10)  ;
%     W = W.* (V * H')  ./ ( W * H * H' + 1e-10) ; 
end




