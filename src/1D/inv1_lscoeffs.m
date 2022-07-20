function coeffs = inv1_lscoeffs(m,L,f0,la,rho)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% L = ifftshit(L);
A = L'*L;
%f0 = 1;

% helpers
fro2 = @(x) norm(x,'fro')^2;
dotM = @(M1,M2) sum( abs(M1'*M2).^2, 'all'); % perform <M1*M1', M2*M2'>
normM2 = @(M) dotM(M,M);
dotT = @(T1,T2) sum( Dnumel(m) .* conj(T1) .* T2, 'all');
normT2 = @(T) sum(Dnumel(m) .* abs(T).^2, 'all');

xx = @(U,TU) 1/la * fro2(L*TU) + fro2(TU) + 1/rho * (normM2(U) - normT2(TU));

xy = @(U,TU,V,TV) 1/la * real((A*TU(:))'*TV(:)) + real(TU(:)'*TV(:)) + ...
	1/rho * (dotM(U,V) - dotT(TU,TV));

%a = @(U,TU,V,TV) 1/2/la * fro2(L*(TU-TV)) + 1/2 * fro2(TU-TV) + ...
%    1/2/rho * (dotM(U,U) + dotM(V,V) - 2*dotM(U,V) - normT2(TU-TV));
%
%b = @(U,TU,V,TV) 1/la * real((A*(TU-TV))'*TV) + real((TU-TV)'*TV) + ...
%    1/rho * (dotM(U,V) - real(dotT(TU,TV)) - dotM(V,V) + normT2(TV));

%c = @(V,TV) 1/2/la*fro2(L*TV) + 1/2*fro2(TV) + 1/2/rho*(dotM(V,V)-normT2(TV));

coeffs = @(U,TU,V,TV) num2cell(1/f0 * [xx(U,TU), xx(V,TV), xy(U,TU,V,TV), 0, 0, 0]);

end
