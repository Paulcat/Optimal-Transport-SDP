function C = Dnumeld(m)
%DNUMELN Number of diagonal elements for Toeplitz matrix

O = ones([m,1]); % 1 only useful for 1d
O = padarray(O,m-1,'post');
C = ifftshift( ifftn(fftn(O).^2) );
