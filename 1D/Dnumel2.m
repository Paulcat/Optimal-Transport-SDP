function C = Dnumel2(n)
%DNUMEL(N) Number of 'diagonal' elements for 2-level Toeplitz matrix

O = ones(n);
O = padarray(O,n-1,'post');
C = ifftshift( ifft2(fft2(O).^2) );

end

