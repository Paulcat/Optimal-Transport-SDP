function [Fn,normalization] = inv_fobj(m,L,la)
%INV_FOBJ Objective for invariant measure problem
%   If T are the moments of the invariant measure (i.e. such that B(m) = m
%   for some operator B), the objective reads
%
%      f(T) = 1/2 * ||T||^2 + 1/(2*la) * ||L*T||^2
%
%   up to some normalization. The matrix L consists in the Fourier
%   coefficients of (B-Id).
%
%   See also INV_FGRAD.



% L = ifftshift(L);
F = @(T) 1/2/la * norm(L*T(:),'fro')^2 + 1/2 * norm(T,'fro')^2;

% normalization
%U0 = ones(m,1)/sqrt(m); %U0(1) = 1; % 1/prod(m)
T0 = ones([2*m-1,1])/prod(m); 

%
normalization = F(T0(:)); % wrt to objective value at U0
%normalization = 1/la*normT(N.*(L'*(L*T0))) + normT(N.*T0) +...
%	1/rho * sqrt(norm(U0'*U0,'fro')^2 - normT2(T0));% wrt to norm of gradient at U0


Fn = @(T) F(T) / normalization;

end
