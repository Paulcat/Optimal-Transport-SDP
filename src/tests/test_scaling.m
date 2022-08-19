% test functional scalings

n = [10,10];
m = n+1;
M = prod(m);

% variable
U = rand(M,8) + 1i*rand(M,8);
v = rand(M,1) + 1i*rand(M,1);

% ot cost
cost = rand(2*m-1);

% moments
c1 = rand(2*n(1)+1,1) + 1i*rand(2*n(1)+1,1);
c2 = rand(2*n(2)+1,1) + 1i*rand(2*n(2)+1,1);

% hyper parameters and scaling constants for ffw
Cl = (norm(eval1(c1),'inf') + norm(eval2(c2),'inf'))/2;
Cr = prod(nn)^2; %??
%
la  = 1e-10*Cl; % "unbalanced" penalization 
rho = 1e-3*Cr; % toeplitz penalization

[f,f0]  = ot1_fobj(m,cost,c1,c2,la,rho);


O = zeros(prod(m),1);

f(O)
f(U)
f(v)