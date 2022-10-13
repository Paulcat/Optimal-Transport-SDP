function fx = fftd(d,x)
%FFTD d-dimensional discrete Fourier transform
%	MYFFT(X) returns the n-dimensional trasnform of (n+1)-dimensional
%	tensor x


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
