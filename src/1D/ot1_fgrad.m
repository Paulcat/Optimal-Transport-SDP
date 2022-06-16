function [G,GU] = ot1_fgrad(m,cost,u1,u2,f0,la,rho)
%FGRAD OT gradient for FFW (1D)
%   G = FGRAD(M,C,U1,U2,LA,RHO) returns the gradient
%
%   [G,P] = FGRAD(...) also returns the gradient with precomputations, for
%   speedup.

debug = 0;

% scaling
%f0 = 4*la / (norm(u1,'fro')^2 + norm(u2,'fro')^2);
% f0 = 1;

% helpers
N = 1./Dnumel2(m);
I1 = zeros(2*m-1); I1(:,1) = 1;
I2 = zeros(2*m-1); I2(1,:) = 1;
I  = I1+I2;
%
Y = u1.*I1 + u2.'.*I2;
%

% gradient part: terms depending on TVALS
GT = @(T,h) 1/2/la * Tprod2(m,N.*(T.*I-Y),h) - 1/rho * Tprod2(m,T,h);

% gradient
G  = @(U,h)   1/f0*( Tprod2(m,N.*cost,h) + GT(Tproj2(m,U),h) + 1/rho*U*(U'*h) );
GU = @(T,U,h) 1/f0*( Tprod2(m,N.*cost,h) + GT(T,h) + 1/rho*U*(U'*h) );

if debug
    disp('debug mode!!');
    GT = @(T,h) 2*Tprod2(m,N.*(T.*I-Y),h);
    G  = @(U,h) Tprod2(m,Tproj2(m,U),h);
end

end

