function [Tvals] = Tproj2(m,U)
%TPROJ2 2-level Toeplitz projection
%   TVALS = TPROJ2(N,U) returns the 'diagonal' values corresponding to the
%   projection of the matrix UU' on the set of 2-level Toeplitz matrices.
%
%   Internal use only

r = size(U,2);
U = reshape(U,[m,r]);

PU = padarray(U,[m-1,0],'post');
Tvals = sum(ifft2(abs(fft2(PU)).^2), 3) ./ Dnumel2(m);

end

