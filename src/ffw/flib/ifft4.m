function ifx = ifft4(x)
%IFFT4 4-dimensional Inverse Fast Fourier Transform
%	Internal use only

ifx = ifft(ifft(ifft2(x),[],3),[],4);

end

