function fx = fftd(d,x)
%FFTD d-dimensional discrete Fourier transform
%   FFTD(X) is the same as FFTN(X) if X is a d-dimensional array. If not,
%   FFTD if the discrete Fourier transform along the first d dimensions.
%
%   Internal use only
%   See also ifftd


if d==1
	fx = fft(x);
elseif d==2
	fx = fft2(x);
elseif d==3
	fx = fft(fft2(x),[],3);
elseif d==4
	fx = fft(fft(fft2(x),[],3),[],4);
else
	error('incorrect dimension of tensor');
end
