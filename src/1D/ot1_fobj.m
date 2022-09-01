function [Fn,normalization] = ot1_fobj(m,cost,u1,u2,la)
%FOBJ OT objective for FFW (1D)
%   Unbalanced optimal transport energy. Marginal constraints are penalized with MMD
%	 If T are the moments of the transport plan, the objective reads
%	 
%	 	f(T) = <cost|T> + 1/(4*la) * (MMD(T1,u1) + MMD(T2,u2))

debug = 0;

% scaling
%f0 = 1/4/la * (norm(u1,'fro')^2 + norm(u2,'fro')^2);
%f0 = 4*la / (norm(u1,'fro')^2 + norm(u2,'fro')^2);
f0 = 1;

% helper
d = length(m);
%normT2 = @(T) sum( Dnumel2(m) .* abs(T).^2, 1:d);

% objective in terms of Toeplitz entries
F = @(T) real(trace(cost'*T)) + ...
	1/4/la * ( norm(T(:,1)-u1,'fro')^2 + norm(T(1,:).'-u2,'fro')^2 );


normalization = f0; %f(0)
%normalization = F(ones(prod(m),1)); % f(1)...

Fn = @(T) F(T) / normalization;

if debug
    warning('debug mode!!');
end

end

