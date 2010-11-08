function [W,H, resid] = NMF(V,k, protos)
% INPUT: V is m*n non-negative  matrix , k is the reduced dimension
% OUTPUT: W is m*k,  H is k*n non-negative matrices  s.t. V ~ WH

% my simple nmf implementation. Don't use it for Factorization anymore as
% the matlab builtin nnmf gives much better result. But still I use it to
% find the loading matrix for given W and V. There the problem with nnmf is
% that also W is changed..

[m,n]    = size(V);
maxiter  = 50;
H        = rand(k,n);

% initialize W with given data for backward usage
if nargin == 3
   if size(protos,1) == m && size(protos,2) ==k
      W = protos;
   else
      error('NMF: W1 does not fit V and k !');
   end
else
   W = rand(m,k);
end


for i= 1:maxiter
   H = H.*(W'*V)./(W' *W * H  + 1e-10)  ;
   if nargin == 2
      W = W.* (V * H')  ./ ( W * H * H' + 1e-10) ;
   end
end

resid = (sum(sum(abs(V - W*H)))/sum(sum(V))) * 100;
