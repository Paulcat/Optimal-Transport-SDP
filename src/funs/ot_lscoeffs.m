function coeffs = ot_lscoeffs(m,cost,u1,u2,f0,la)
%OT_LSCOEFFS Linesearch coefficients for mmd-penalized ot
%   OT_LSCOEFFS(...) returns the coefficients of the quadratic form
%
%      q(a,b) = f(a*T + b*S) = xx*a^2 + xy*a*b + yy*b^2 + x*a + y*b + cst
%
%   where f is the ot objective with penalized marginal constraints.
%
%   See also OT_FOBJ.

debug = 0;

d = length(m); % dimension of the plan
dm = d/2; % dimension of marginals

% subroutines
fro2 = @(x) norm(x(:),'fro')^2;
dotp = @(x,y) real(x(:)'*y(:));
%
I = cell(1,d); 
I(1:dm) = {':'}; I(dm+1:d) = {1};
J = flip(I);
%
C1 = @(T) T(I{:}); % T(:,:,...,1,1,...) in arbitrary dimension
C2 = @(T) T(J{:});



% coefficients
xx = @(U,TU) 1/2/la * (fro2(C1(TU)) + fro2(C2(TU))) ;
xy = @(U,TU,V,TV) 1/2/la * (dotp(C1(TU),C1(TV)) + dotp(C2(TU),C2(TV)));
x 	= @(T) dotp(cost,T) - 1/2/la * (dotp(C1(T),u1) + dotp(C2(T),u2));
z  = 1/4/la * (fro2(u1) + fro2(u2));
%
coeffs = @(U,TU,V,TV) num2cell(1/f0 * [xx(U,TU),xx(V,TV),xy(U,TU,V,TV),x(TU),x(TV),z]); %TODO: why 1/f0 at the end??



if debug
	disp('debug mode!!');
	xx = @(U,TU) 1/2/la * (fro2(C1(TU)));
	xy = @(U,TU,V,TV) 1/2/la*real(C1(TU)'*C1(TV));
	x  = @(T) -1/2/la*real(C1(T)'*u1(:));
	coeffs = @(U,TU,V,TV) num2cell( [xx(U,TU),xx(V,TV),xy(U,TU,V,TV),x(TU),x(TV),1/4/la*norm(u1,'fro')^2] );
end

end

