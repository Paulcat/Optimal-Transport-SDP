function [mu,nu,stop] = ffw_ls(ls_co,U,T,v,s,f0,rho,Dnumel,flag)
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here

debug = 0;
stop = 0;


dotT = @(T1,T2) real(sum(Dnumel .* conj(T1) .* T2, 'all'));
normT = @(T) sum(Dnumel .* abs(T).^2, 'all');
dotM = @(M1,M2) real(sum(abs(M1'*M2).^2, 'all'));

co = ls_co(U,T,v,s);
[cxx,cyy,cxy,cx,cy,cst] = deal(co{:});
%
cxx = cxx + 1/rho/f0 * (dotM(U,U) - normT(T));
cxy = cxy + 1/rho/f0 * (dotM(U,v) - dotT(T,s));
cyy = cyy + 1/rho/f0 * (dotM(v,v) - normT(s));

%TODO: unbounded linesearch???

if flag == 1 % linear constraint is mu+nu=1, e.g. for invariant
	a = cxx+cyy-2*cxy;
	if a == 0 % can only happen at first iteration. Then S=X and we have found an optimum?
		mu = 1;
		nu = 0;
		stop = 1;
		
	else % mu + nu = 1 and S != X
		b = cxy-cyy+cx-cy; % we solve min a/2*mu^2 + b*mu
		%mu = min( max(0,-cx/(2*cxx)), 1);
		mu = min( max(0, -b/a), 1);
		nu = 1 - mu;
	end

	return;
end


%denom = 4*cxx*cyy - cxy^2;
%num1  = cxy*cy - 2*cyy*cx;
%num2  = cxy*cx - 2*cxx*cy;
det   = cxx*cyy - cxy^2;
num_a = cxy*cy - cyy*cx; %or opposite??
num_b = cxy*cx - cxx*cy; %or opposite?

if det == 0 % X=0 or S=tX, hence linesearch along rays
    if cxx == 0 % X=0, look along ray S
		 mu = 0;
		 nu = max(0,-cy/cyy); % unbounded linesearch
	 else % S=tX, look along ray X
	 	 mu = max(0,-cx/cxx);
		 nu = 0;
    end
else
	mu = num_a/det;
	nu = num_b/det;
    
	%if ~(mu>=0) || ~(nu>=0)
	%	f1 = cyy*nu^2 + cy*nu; % mu = 0?
	%	f2 = cxx*mu^2 + cx*mu; % or nu = 0?
	%
	%	I = f1 < f2; % if I then mu=0, else nu=0
	%	mu = (1-I)*max(0,mu);
	%	nu = I*max(0,nu);
   %end
	if (mu >= 0) && (nu < 0)
		mu = max(0,-cx/cxx);
		nu = 0;
	elseif (mu < 0)
		mu = 0;
		nu = max(0,-cy/cyy);
	end
end

if debug
    A = 1/2 * [cxx, cxy; cxy, cyy];
    B = [cx, cy];
    warning off
    S = sqrtm(A);
    %warning on
    %cvx_solver mosek
    cvx_precision high
    cvx_begin %quiet
    variable X(2)
    %
    X >= 0
    %sum(X) <= 1
    %
    %minimize( norm(S*X,'fro')^2 + B*X )
	 minimize( X'*A*X + B*X + cst)
    cvx_end
    
    err1 = abs(X(1)-mu)/abs(X(1));
    err2 = abs(X(2)-nu)/abs(X(2));
    fprintf('\t linesearch error 1: %d\n', err1);
    fprintf('\t linesearch error 2: %d\n', err2);
end

end

