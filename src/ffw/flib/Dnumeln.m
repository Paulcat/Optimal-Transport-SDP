function C = Dnumeln(m)
%DNUMELN Number of diagonal elements for n-level Toeplitz matrix
%   DNUMELN(M) returns a tensor D such that D(k1,...,kn) contains the
%   number of tuples (i,j) satisfying
%
%      T(i+k,l) = T(i,l-k)
%
%   in a n-level Toeplitz matrix of size (M1 x ... x Mn).
%
%   For instance, if n=1, D is a vector containing the number of elements
%   on each diagonal.

O = ones([m,1]); % 1 only useful for 1d
O = padarray(O,m-1,'post');
C = ifftshift( ifftn(fftn(O).^2) );
