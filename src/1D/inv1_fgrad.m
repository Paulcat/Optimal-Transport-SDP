function [G,GU] = inv1_fgrad(m,L,f0,la,rho)
%INV1_FGRAD

% A = ifftshift(L'*L);
A = L'*L;
N = 1./Dnumel(m);

GT = @(T,h) 1/la*Tprod(m,N.*(A*T),h) + Tprod(m,N.*T,h) - 1/rho*Tprod(m,T,h);

G  = @(U,h) 1/f0 * ( 1/rho * U*(U'*h) + GT(Tproj(m,U),h) );
GU = @(T,U,h) 1/f0 * ( 1/rho * U*(U'*h) + GT(T,h) );
% GU = @(T,U,h) 1/la*Tprod(m,N.*(A*T),h) + Tprod(m,N.*T,h);

end