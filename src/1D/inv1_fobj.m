function [Fn,normalization] = inv1_fobj(m,L,la)
% INV1_FOBJ objective for invariant measure problem

d = length(m);
if d ~= 1
	error('dimension should be 1');
end
% L = ifftshift(L);

%normT2 = @(T) sum(Dnumel(m) .* abs(T).^2,1:d);
%normT  = @(T) sqrt(normT2(T));

%FT = @(T) 1/2/la * norm(L*T,'fro')^2 + 1/2 * norm(T,'fro')^2 - 1/2/rho * normT2(T);
%F  = @(U) FT(Tproj(m,U)) + 1/2/rho * norm(U'*U,'fro')^2;


F = @(T) 1/2/la * norm(L*T,'fro')^2 + 1/2 * norm(T,'fro')^2;

% normalization
%U0 = ones(m,1)/sqrt(m); %U0(1) = 1; % 1/prod(m)
T0 = ones(2*m-1,1)/m; 
%N  = 1./Dnumel(m);
%
normalization = F(T0); % wrt to objective value at U0
%normalization = 1/la*normT(N.*(L'*(L*T0))) + normT(N.*T0) +...
%	1/rho * sqrt(norm(U0'*U0,'fro')^2 - normT2(T0));% wrt to norm of gradient at U0


Fn = @(T) F(T) / normalization;

end
