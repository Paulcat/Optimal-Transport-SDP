function Tvals = Tprojn(m,U)
%TPROJn n-level Toeplitz projection
%	TVALS = TPROJN(m,U) returns the 'diagonal' values corresponding to the
%	projection of the matrix UU* on the set of n-level Toeplitz matrices.
%
%	TVALS is a nD-array. Colexicographical ordering is assumed, and the 
%	level-0 moment is the first entry in TVALS (fft-style),
%
%	Internal use only

d = length(m);
r = size(U,2);
U = reshape(U,[m,r]);

PU = padarray(U,[m-1,0],'post');
Tvals = sum( ifftd(d,abs(fftd(d,PU)).^2), 3) ./ Dnumeln(m);

end
