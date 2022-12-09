function coeffs = inv_lscoeffs(m,L,f0,la)
%INV_LSCOEFFS Linesearch coeffcients for invariant problem
%   INV_LSCOEFFS(...) returns the coefficients of the quadratic form
%
%      q(a,b) = f(a*T + b*S) = 1/2*xx*a^2 + xy*a*b + 1/2*yy*b^2 + x*a + y*b + cst
%
%   where f is the invariant problem objective.
%
%   See also INV_FOBJ.




% L = ifftshit(L);
%A = L'*L;
%f0 = 1;

% helpers
fro2 = @(x) norm(x,'fro')^2;
dotp = @(x,y) real(x(:)'*y(:));


% coefficients
xx = @(U,TU) 1/la * fro2(L*TU(:)) + fro2(TU);
%xy = @(U,TU,V,TV) 1/la * real((A*TU(:))'*TV(:)) + real(TU(:)'*TV(:)) ;
xy = @(U,TU,V,TV) 1/la * dotp(L*TU(:),L*TV(:)) + dotp(TU,TV);
%
coeffs = @(U,TU,V,TV) num2cell(1/f0 * [xx(U,TU), xx(V,TV), xy(U,TU,V,TV), 0, 0, 0]);

end
