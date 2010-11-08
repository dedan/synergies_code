function test()

x = [1 2 3; 1 2 3; 8 8 7];
null(x)
orthcomp(x')
orthcomp(null(x))
orthcomp(orthcomp(x'))



function BCOMP = orthcomp(B)

% orthcomp  Orthogonal complement of a subspace.
%
% BCOMP = orthcomp(B) returns a basis for the 
% orthogonal complement of the column space of B.
% This subspace contains all vectors orthogonal
% to the column space of B.
% It is the left nullspace of B.
%
% See also leftnull, nulbasis.

BCOMP = leftnull(B);


function LN = leftnull(A)

% leftnull  Basis for the left nullspace.
%
% LN = leftnull(A) returns a basis for the 
% left nullspace in the *columns* of LN.
%
% The left nullspace of A is the nullspace of A'.
% The command fourbase(A) finds a different basis
% for the left nullspace of A. 
%
% See also fourbase.

LN = nulbasis(A');



function N = nulbasis(A)

% nulbasis  Basis for nullspace.
%
% N = nulbasis(A) returns a basis for the nullspace of A
% in the columns of N. The basis contains the n-r special 
% solutions to Ax=0.  freecol is the list of free columns.
%
% Example:
%
% >> A = [1 2 0 3;
%        [0 0 1 4];
%
% >> N = nulbasis(A)
%
%    N = [-2  -3]   
%        [ 1   0]
%        [ 0  -4]
%        [ 0   1]
%
% See also fourbase.

[R, pivcol] = rref(A, sqrt(eps));
[m, n] = size(A);
r = length(pivcol);
freecol = 1:n;
freecol(pivcol) = [];
N = zeros(n, n-r);
N(freecol, : ) = eye(n-r);
N(pivcol,  : ) = -R(1:r, freecol);
