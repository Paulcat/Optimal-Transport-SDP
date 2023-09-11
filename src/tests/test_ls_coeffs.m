% simple test to verify linesearch coefficients are correct. Also useful to
% check which kind of inputs our functionals get.

problem = 'invariant';

n = [10,10];
m = n+1;
M = prod(m);

% variable
U = rand(M,8) + 1i*rand(M,8);
v = rand(M,1) + 1i*rand(M,1);
%
TU = Tprojn(m,U);
Tv = Tprojn(m,v);

% hyperparameters
la  = rand;

switch problem
	case 'OT'
		% moments
		c1 = rand(2*n(1)+1,1) + 1i*rand(2*n(1)+1,1);
		c2 = rand(2*n(2)+1,1) + 1i*rand(2*n(2)+1,1);
		
		% ot cost
		cost = rand(2*m-1);
		
		% functionals
		[f,f0]  = ot_fobj(m,cost,c1,c2,la);
		ls		  = ot_lscoeffs(m,cost,c1,c2,f0,la);
		
	case 'invariant'
		L = rand(prod(2*n+1)) + 1i * rand(prod(2*n+1));
		
		% functionals
		[f,f0] = inv_fobj(m,L,la);
		ls		 = inv_lscoeffs(m,L,f0,la);
end

co		  = ls(U,TU,v,Tv);


[cxx,cyy,cxy,cx,cy,cst] = deal(co{:});
A = 1/2*[cxx,cxy;cxy,cyy];
B = [cx,cy];
C = cst;

[mu,nu] = deal(rand,rand);
al = [mu;nu];
W = [sqrt(mu)*U, sqrt(nu)*v];
T = Tprojn(m,W);


test1 = f(T);
test2 = al'*A*al + B*al + C;
fprintf('mismatch: %d\n', abs(test1-test2)/abs(test1));


