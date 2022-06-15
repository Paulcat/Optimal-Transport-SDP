% simple test to verify linesaerch coefficients are correct. Also useful to
% check which kind of inputs our functionals get.

n = [10,10];
m = n+1;
M = prod(m);

% variable
U = rand(M,8) + 1i*rand(M,8);
v = rand(M,1) + 1i*rand(M,1);

% moments
c1 = rand(2*n(1)+1,1) + 1i*rand(2*n(1)+1,1);
c2 = rand(2*n(2)+1,1) + 1i*rand(2*n(2)+1,1);

% hyperparameters
la  = rand;
rho = rand;

% ot cost
cost = rand(2*m-1);

% functionals
[f,f0]  = ot1_fobj(m,cost,c1,c2,la,rho);
p		  = ot1_lscoeffs(m,cost,c1,c2,la,rho);
co		  = p(U,Tproj2(m,U),v,Tproj2(m,v));

[cxx,cyy,cxy,cx,cy,cst] = deal(co{:});
A = 1/2 * [cxx,cxy;cxy,cyy];
B = [cx,cy];
C = cst;

[mu,nu] = deal(rand,rand);
t = [mu;nu];
W = [sqrt(mu)*U, sqrt(nu)*v];


test1 = f(W);
test2 = t'*A*t + B*t + C;
fprintf('mismatch: %d\n', abs(test1-test2)/abs(test1));


