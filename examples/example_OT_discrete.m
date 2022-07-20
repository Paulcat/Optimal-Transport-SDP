% example OT between 1D discrete measures


% ***** set up problem ******
% ***************************
% create marginals
[s1,s2]  = deal(4,3);
% 1st
x1 = rand(s1,1);
a1 = ones(s1,1)/s1;
% 2nd
x2 = rand(s2,1);
a2 = ones(s2,1)/s2;


% moments of marginals
[n1,n2]  = deal(10,10);
nn = [n1 n2];
% 1st marginal
c1 = exp(-2i*pi*(-n1:n1)'*x1(:)') * a1(:);
F1 = exp(-2i*pi*(0:n1)'*x1(:)');
U1 = F1 .* sqrt(a1');
% 2nd marginal
c2 = exp(-2i*pi*(-n2:n2)'*x2(:)') * a2(:);
F2 = exp(-2i*pi*(0:n2)'*x2(:)');
U2 = F2 .* sqrt(a2');
% shift (order 0 comes first)
c1 = ifftshift(c1);
c2 = ifftshift(c2);


% polynomial display
[N1,N2] = deal(512,512);
[t1,t2] = deal( (0:N1-1)'/N1, (0:N2-1)'/N2 );
NN 	  = [N1,N2];
%
eval1 = @(p) real( exp(-2i*pi*(-n1:n1)'*t1')' * p(:) );
eval2 = @(p) real( exp(-2i*pi*(-n2:n2)'*t2')' * p(:) );


% hyper parameters and scaling constants for ffw
Cl = (norm(eval1(c1),'inf') + norm(eval2(c2),'inf'))/2;
%Cr = prod(nn)^2; %??
Cr = prod(nn) * (norm(c1,'inf') + norm(c2,'inf'))/2;
%
la  = 1e-8*Cl; % "unbalanced" penalization 
rho = 1e-3*Cr; % toeplitz penalization



% transport cost (in Fourier domain)
cost = zeros(2*nn+1);
cost(n1:n1+2,n2:n2+2) = [0 0 -1/4; 0 1/2 0; -1/4 0 0];
% shift
cost = ifftshift(cost);


% variable's size and functionals
mm		 = nn+1;
[f,f0] = ot1_fobj(mm,cost,c1,c2,la,rho);
[g,gU] = ot1_fgrad(mm,cost,c1,c2,f0,la,rho);

% load problem
problem = struct;
%
problem.name		= 'OT';
problem.vardim 	= mm; % size
problem.fobj 		= f; % objective
problem.gscaling	= 1/f0; % scaling constant st f(0)=1
problem.grad 		= g; % gradient
problem.grad_pre  = gU; % gradient with partial precomputations
problem.hparams	= [la,rho];
problem.ls 			= ot1_lscoeffs(mm,cost,c1,c2,f0,la,rho); % coefficients for line-search
% ******************
% ******************




% ***** load solver options *****
% ********************************
options = struct;
%
%
options.tol 			= 1e-5; % tolerance on ffw criterion
options.maxiter		= 6; % max iterations for ffw
options.bfgsProgTol 	= 1e-16; % tolerance on ?
options.bfgsMaxIter 	= 500; 
options.lmoTol 		= 1e-10;
options.lmoMaxIter 	= 1e3;
% ************************
% ************************



% *** cvx ground truth ***
% ************************
% periodic cost
C = sin(x1(:)-x2(:)').^2;

% cvx solver
cvx_begin
variable P0(s1,s2) nonnegative
%
sum(P0,1) == a2.';
sum(P0,2) == a1;
%
minimize( trace(C(:)'*P0(:)) );
cvx_end

% with Hungarian algorithm
%rho = hungarianLSAP(C);

% recovery
[I,J] = find(P0 > 1e-7);




% ********************
% ********************


% *** OUR SOLVER ***
% ******************
U = FFW(problem,options);

% prony extraction
options_prony.factorized = 1;
[x,a] = mvprony(U,nn,2,options_prony);

% display
clf, hold on;
% 'ground truth'
scatter(x1(I),x2(J),'o');
% ffw + prony
scatter(x(:,1),x(:,2),50,'x','linewidth',2);
scatter(x1,zeros(s1,1),100,'.');
scatter(zeros(s2,1),x2,100,'.');
xlim([0,1]), ylim([0,1]);
% ********************
% ********************
