function [Fn,normalization] = inv2_fobj(m,L,la,rho)
% INV"_FOBJ objective for invariant measure problem

d = length(m);
if d ~= 2
	error('dimension should be 2');
end

normT2 = @(T) sum(Dnumel(m) .* abs(T).^2, 1:d);

FT = @(T) 1/2/la * norm(L*T,'fro')^2 + 1/2 * norm(T,'fro') - 1/2/rho * normT2(T);
F  = (U) FT(Tproj2(m,U)) + 1/2/rho * norm(U'*U,'fro')^2;

% normalization
U0 = ones(prod(m),1)/sqrt(prod(m));
normalization = F(U0);

Fn = @(U) F(U)/normalization;

end
