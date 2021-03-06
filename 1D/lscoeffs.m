function coeffs = lscoeffs(m,cost,u1,u2,la,rho)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here

% scaling
f0 = 4*la / (norm(u1,'fro')^2 + norm(u2,'fro')^2);

% helpers
fro2   = @(x) norm(x,'fro')^2;
normT2 = @(T) sum( Dnumel2(m) .* abs(T).^2, 'all');
dotT   = @(T1,T2) sum( Dnumel2(m) .* conj(T1) .* T2, 'all');
dotM   = @(M1,M2) sum( abs(M1'*M2).^2, 'all');
C1     = @(T) T(:,1);
C2     = @(T) T(1,:).';

%
a = @(U,TU) 1/4/la * (fro2(C1(TU)) + fro2(C2(TU))) + ...
    1/2/rho * (dotM(U,U) - normT2(TU));
%
b = @(U,TU,V,TV) 1/2/la * (real(C1(TU)'*C1(TV) + C2(TU)'*C2(TV))) + ...
    1/rho * (dotM(U,V) - real(dotT(TU,TV)));
%
c = @(T) real( trace(cost'*T) - 1/2/la * (C1(T)'*u1(:) + C2(T)'*u2(:)) );
%

coeffs = @(U,TU,V,TV) num2cell(f0 * [a(U,TU),a(V,TV),b(U,TU,V,TV),c(TU),c(TV)]);
end

