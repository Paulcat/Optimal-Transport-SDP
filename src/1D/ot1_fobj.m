function [Fn,normalization] = ot1_fobj(m,cost,u1,u2,la,rho)
%FOBJ OT objective for FFW (1D)
%   Detailed explanation goes here

debug = 0;

% scaling
f0 = 1/4/la * (norm(u1,'fro')^2 + norm(u2,'fro')^2);
%f0 = 4*la / (norm(u1,'fro')^2 + norm(u2,'fro')^2);
% f0 = 1;

% helper
d = length(m);
normT2 = @(T) sum( Dnumel2(m) .* abs(T).^2, 1:d);

% objective part: terms depending on TVALS
FT = @(T) real(trace(cost'*T)) + ...
    1/4/la  * ( norm(T(:,1)-u1,'fro')^2 + norm(T(1,:).'-u2,'fro')^2 ) + ...
    -1/2/rho * normT2(T);

% objective
%F = @(U) f0 * ( FT(Tproj2(m,U)) + 1/2/rho * norm(U'*U,'fro')^2 );
%normalization = f0;
F = @(U) FT(Tproj2(m,U)) + 1/2/rho * norm(U'*U, 'fro')^2;

normalization = f0; %f(0)
%normalization = F(ones(prod(m),1)); % f(1)...

Fn = @(U) F(U) / normalization;

if debug
    warning('debug mode!!');
    FT = @(T) 1/4/la * ( norm(T(:,1)-u1,'fro')^2 );
    %F = @(U) 1/2*normT2(Tproj2(m,U));
	 F = @(U) FT(Tproj2(m,U));
end

end

