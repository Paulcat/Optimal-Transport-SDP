function [Fn,normalization] = ot_fobj(m,cost,u1,u2,la)
%OT_FOBJ Objective for mmd-penalized ot (in arbitrary dimension)
%   Unbalanced optimal transport energy. Marginal constraints are penalized
%   with MM discrepancies: if T are the moments of the transport plan, the
%   objective reads
%	 
%	 	f(T) = <cost,T> + 1/(4*la) * (MMD(T1,u1) + MMD(T2,u2))
%
%	 up to some normalization.
%
%	 See also OT_FGRAD.

debug = 0;


% scaling
%f0 = 1/4/la * (norm(u1,'fro')^2 + norm(u2,'fro')^2);
f0 = 1;


% helpers
flat = @(x) x(:);
%
d 	= length(m); % dimension of transport plan
dm = d/2; % dimension of marginals
%
i1 = cell(1,d);
i1(1:dm) = {':'}; i1(dm+1:d) = {1}; % indices corresponding to 1st marginal
i2 = flip(i1); % to 2nd marginal


% normalized objective
F = @(T) (	real(cost(:)'*T(:)) + ... % OT cost
	1/4/la * ( 	norm(flat(T(i1{:}))-u1(:),'fro')^2 + ... % marginals penalization
					norm(flat(T(i2{:}))-u2(:),'fro')^2 ) );

normalization = f0;
%normalization = F(ones(prod(m),1)); % f(1)...
Fn = @(T) F(T) / normalization;



if debug
    warning('debug mode!!');
end

end

