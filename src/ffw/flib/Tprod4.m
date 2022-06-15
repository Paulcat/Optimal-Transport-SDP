function Tx = Tprod4(m,Tvals,x)
%TPROD4 4-level Toeplitz product
%   TPROD4(M,TVALS,X) returns the product between the 4-level Toeplitz
%   matrix whose 'diagonal' values are TVALS, and the array X.
%
%   M = (M1,M2,M3,M4) is the dimension of the 4-level Toeplitz matrices.
%
%   TVALS is a tensor of size (2*M1-1,2*M2-1,2*M3-1,2*M4-1)
%   TODO: explain order
%
%   X is a matrix of size (M1*M2*M3*M4,r).
%
%   Internal use only.

r = size(x,2);
x = reshape(x,[m,r]);

Px = padarray(x,[m-1,0],'post');
Cx = ifft4( fft4(Tvals) .* fft4(Px) ); % broadcasting...

Tx = Cx(1:m(1),1:m(2),1:m(3),1:m(4),:);
Tx = reshape(Tx,prod(m),r);

end

