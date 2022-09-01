function coeffs = ot1_lscoeffs(m,cost,u1,u2,f0,la,rho)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here

debug = 0;


% scaling
%f0 = 4*la / (norm(u1,'fro')^2 + norm(u2,'fro')^2);

% helpers
fro2   = @(x) norm(x,'fro')^2;
%normT2 = @(T) sum( Dnumel2(m) .* abs(T).^2, 'all');
%dotT   = @(T1,T2) sum( Dnumel2(m) .* conj(T1) .* T2, 'all');
%dotM   = @(M1,M2) sum( abs(M1'*M2).^2, 'all');
C1     = @(T) T(:,1);
C2     = @(T) T(1,:).';

%
xx = @(U,TU) 1/2/la * (fro2(C1(TU)) + fro2(C2(TU))) ;%+ ... %before it was 1/4/la...
%	1/rho * (dotM(U,U) - normT2(TU));
%
xy = @(U,TU,V,TV) 1/2/la * (real(C1(TU)'*C1(TV) + C2(TU)'*C2(TV))) ;%+ ...
%	1/rho * (dotM(U,V) - real(dotT(TU,TV)));
%
x = @(T) real( trace(cost'*T) - 1/2/la * (C1(T)'*u1(:) + C2(T)'*u2(:)) ); % before it was 1/2/la...
%

a = 1/f0;

coeffs = @(U,TU,V,TV) num2cell(1/f0 * [xx(U,TU),xx(V,TV),xy(U,TU,V,TV),x(TU),x(TV),a]);

if debug
	disp('debug mode!!');
	xx = @(U,TU) 1/2/la * (fro2(C1(TU)));
	xy = @(U,TU,V,TV) 1/2/la*real(C1(TU)'*C1(TV));
	x  = @(T) -1/2/la*real(C1(T)'*u1(:));
	coeffs = @(U,TU,V,TV) num2cell( [xx(U,TU),xx(V,TV),xy(U,TU,V,TV),x(TU),x(TV),1/4/la*norm(u1,'fro')^2] );
end

end

