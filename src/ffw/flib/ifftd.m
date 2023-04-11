function ifx = ifftd(d,x)
%IFFTD d-dimensional inverse discrete Fourier transform
%   IFFTD(D,X) is the same as IFFTN(X) if X is a d-dimensional array. If
%   not, IFFTD performs the discrete inverse Fourier transform along the
%   first d dimensions.
%
%   See also fftd

if d==1
	ifx = ifft(x);
elseif d==2
	ifx = ifft2(x);
elseif d==3
	ifx = ifft(ifft2(x),[],3);
elseif d==4
	ifx = ifft(ifft(ifft2(x),[],3),[],4);
else
	error('dimension of tensor is incorrect');
end
