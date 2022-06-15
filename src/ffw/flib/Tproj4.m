function [Tvals] = Tproj4(m,U)
%TPROJ4 4-level Toeplitz projection
%   TVALS = TPROJ4(M,U) returns the 'diagonal' values corresponding to the
%   projection of the matrix UU' on the set of 4-level Toeplitz matrices.
%
%   TVALS is tensor-shaped. Colexicographical ordering is assumed, and the
%   level-0 moment is the first entry in TVALS (fft-style).
%
%   Internal use only

r = size(U,2);
U = reshape(U,[m,r]);

PU = padarray(U,[m-1,0],'post');
Tvals = sum(ifft4(abs(fft4(PU)).^2), 5) ./ Dnumel4(m);

end

