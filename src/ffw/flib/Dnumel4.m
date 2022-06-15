function C = Dnumel4(m)
%UNTITLED11 Summary of this function goes here
%   Detailed explanation goes here

O = ones(m);
O = padarray(O,m-1,'post');
C = ifftshift( ifftn(fftn(O).^2) );

end

