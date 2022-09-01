function coeffs = inv2_lscoeffs(m,L,f0,la)
%UNTITLED

A = L'*L;

fro2 	 = @(x) 		norm(x,'fro')^2;
%dotM 	 = @(M1,M2) sum( abs(M1'*M2).^2, 'all');
%normM2 = @(M) 		dotM(M,M);
%dotT 	 = @(T1,T2) sum( Dnumel2(m) .* conj(T1) .* T2, 'all');
%normT2 = @(T) 		sum( Dnumel2(m) .* abs(T).^2, 'all');

xx = @(U,TU) 1/la * fro2(L*TU) + fro2(TU) ;%+ 1/rho * (normM2(U) - normT2(TU));

xy = @(U,TU,V,TV) 1/la * real((A*TU(:))'*TV(:)) + real(TU(:)'*TV(:)) ;%+ ...
	%1/rho * (dotM(U,V) - dotT(TU,TV));

coeffs = @(U,TU,V,TV) num2cell(1/f0 * [xx(U,TU), xx(V,TV), xy(U,TU,V,TV), 0, 0, 0]);

end
