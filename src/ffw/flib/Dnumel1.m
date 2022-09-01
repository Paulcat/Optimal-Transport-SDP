function C = Dnumel(m)
%DNUMEL(N) Number of diagonal elements for Toeplitz matrix

O = ones(m,1);
O= padarray(O,m-1,'post');
C = ifftshift( ifft(fft(O).^2) );

end
