function [G,GU] = inv2_fgrad(m,L,f0,la,rho)
%INV2_FGRAD

A = L'*L;
N = 1./Dnumel2(m);

GT = @(T,h) 	1/la*Tprod(m,N.*(A*T),h) + Tprod2(m,N.*T,h) - 1/rho*Tprod2(m,T,h);
G  = @(U,h) 	1/f0 * ( 1/rho * U*(U'*h) + GT(Tproj2(m,U),h) );
GU = @(T,U,h)	1/f0 * ( 1/rho * U*(U'*h) + GT(T,h) );

end

