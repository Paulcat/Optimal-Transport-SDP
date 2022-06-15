function [fx] = fft4(x)
%FFT4 4-dimensional Fast Fourier Transform
%	Internal use only

fx = fft(fft(fft2(x),[],3),[],4);

end

