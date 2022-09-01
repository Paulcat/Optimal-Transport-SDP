function [G] = inv1_fgrad(m,L,f0,la)
%INV1_FGRAD

% A = ifftshift(L'*L);
A = L'*L;
N = 1./Dnumel1(m);

%GT = @(T,h) 1/la*Tprod(m,N.*(A*T),h) + Tprod(m,N.*T,h) - 1/rho*Tprod(m,T,h);

%G  = @(U,h) 1/f0 * ( 1/rho * U*(U'*h) + GT(Tproj(m,U),h) );
%GU = @(T,U,h) 1/f0 * ( 1/rho * U*(U'*h) + GT(T,h) );
% GU = @(T,U,h) 1/la*Tprod(m,N.*(A*T),h) + Tprod(m,N.*T,h);

G = @(T,h) 1/f0 * Tprod1(m,N.*( 1/la*(A*T) + T),h);

end
