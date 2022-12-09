function Tx = Tprodn(m,Tvals,x)
%TPRODN n-level Toeplitz product
%   TPRODN(M,TVALS,X) returns the product between the n-level Toeplitz
%   matrix specified by TVALS and the array X.
%
%   M = (M1,...,Mn) gives the dimensions of the Toeplitz matrix.
%
%   TVALS is a tensor of size (2*M1-2,...,2*Mn-1) containing all 'diagonal'
%   values, in fft order (ie first entry is the main diagonal = 0-level
%   moment).
%
%   X is a matrix of size (M1 x ... x Mn x r).
%
%   Internal use only.
%   See also Tprojn.


d = length(m);
r = size(x,2);
x = reshape(x,[m,r]);

Px = padarray(x,[m-1,0],'post');
Cx = ifftd(d, fftd(d,Tvals) .* fftd(d,Px) ); % broadcasting

% cropping in arbitrary dimension requires cells...
crop = cellfun(@(j) 1:j, num2cell(m), 'un', 0);
Tx = Cx(crop{:},:);
Tx = reshape(Tx,prod(m),r);
