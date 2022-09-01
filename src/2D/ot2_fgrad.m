function [G,GU] = ot2_fgrad(m,cost,u1,u2,f0,la)
%FGRAD OT gradient for FFW (2D)

debug = 0;

% scaling
%f0 = 4*la / (norm(u1(:),'fro')^2 + norm(u2(:),'fro')^2);
%f0 = 1;

% helpers
N = 1./Dnumel4(m);
I1 = zeros(2*m-1); I1(:,:,1,1) = 1;
I2 = zeros(2*m-1); I2(1,1,:,:) = 1;
I  = I1 + I2;
%
Y  = u1.*I1 + reshape(u2,[1,1,size(u2)]).*I2;
%

% gradient part: terms depending on TVALS
%GT = @(T,h) 1/2/la * Tprod4(m,N.*(T.*I-Y),h);

% gradient
%G  = @(U,h)   f0*( Tprod4(m,N.*cost,h) + GT(Tproj4(m,U),h));
%GU = @(T,U,h) f0*( Tprod4(m,N.*cost,h) + GT(T,h));
G = @(T,h) 1/f0 * ( Tprod4(m,N.*cost,h) + 1/2/la * Tprod4(m,N.*(T.*I-Y),h) );

if debug
    disp('debug mode!!');
    G = @(U,h) Tprod4(m,Tproj4(m,U),h);
end

end

