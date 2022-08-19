% example of FFW algorithm for OT between gaussian measures


% ***** set up problem ******
% ***************************
% create marginals
s = 6;
offset = 0.01;
% 1st
x1 = linspace(0,1,s)';
a1 = offset + exp(-(x1-.3).^2/2/.01^2);
a1 = a1/norm(a1,1);
% 2nd
x2 = linspace(0,1,s)';
a2 = offset + .3*exp(-(x2-.7).^2/2/.005^2) + .6*exp(-(x2-.5).^2/2/.02^2);
a2 = a2/norm(a2,1);

% display
clf, hold on;
area(x1,a1);
area(x2,a2);


% moments of marginals
[n1,n2] = deal(20,20);
nn 	  = [n1,n2];
% 1st marginal
c1 = exp(-2i*pi*(-n1:n1)'*x1(:)') * a1(:);
F1 = exp(-2i*pi*(0:n1)'*x1(:)');
U1 = F1 .* sqrt(a1');
% 2nd marginal
c2 = exp(-2i*pi*(-n2:n2)'*x2(:)') * a2(:);
F2 = exp(-2i*pi*(0:n2)'*x2(:)');
U2 = F2 .* sqrt(a2');
% shift
c1 = ifftshift(c1);
c2 = ifftshift(c2);

% polynomial display
[N1,N2] = deal(512,512);
[t1,t2] = deal( (0:N1-1)'/N1, (0:N2-1)'/N2 );
NN      = [N1,N2];
%
eval1 = @(p) real( exp(-2i*pi*(-n1:n1)'*t1')' * p(:) );
eval2 = @(p) real( exp(-2i*pi*(-n2:n2)'*t2')' * p(:) );

% display
clf, hold on;
plot(t1,eval1(fftshift(c1)));
plot(t2,eval2(fftshift(c2)));


% hyper parameters and scaling constants for ffw
Cl = (norm(eval1(c1),'inf') + norm(eval2(c2),'inf'))/2;
Cr = prod(nn)^2; %??
%
la  = 1e-3*Cl; % "unbalanced" penalization
rho = 1e-2*Cr; % toeplitz penalization


% transport cost (in Fourier domain)
cost = zeros(2*nn+1);
cost(n1:n1+2,n2:n2+2) = [0 0 -1/4; 0 1/2 0; -1/4 0 0];
% shift
cost = ifftshift(cost);



% variable's size and functionals
mm     = nn+1;
[f,f0] = ot1_fobj(mm,cost,c1,c2,la,rho);
[g,gU] = ot1_fgrad(mm,cost,c1,c2,f0,la,rho);

% load problem
problem = struct;
%
problem.name 		= 'OT';
problem.vardim 	= mm;
problem.fobj 		= f;
problem.gscaling 	= 1/f0;
problem.grad 		= g;
problem.grad_pre 	= gU;
problem.hparams 	= [la,rho];
problem.ls 			= ot1_lscoeffs(mm,cost,c1,c2,f0,la,rho);
% ***************
% ***************

% **** load solver options *****
% ******************************
options = struct;
%
options.tol 			= 1e-5;
options.maxiter 		= 35;
options.bfgsProgTol 	= 1e-16;
options.bfgsMaxIter 	= 500;
options.lmoTol 		= 1e-10;
options.lmoMaxIter 	= 1e3;
% ***********************
% ***********************



% **** cvx ground truth ****
% **************************
% periodic cost
C = sin(x1(:)-x2(:)').^2;

% cvx solver
cvx_begin
variable P0(s,s) nonnegative
%
sum(P0,1) == a2.';
sum(P0,2) == a1;
%
minimize( trace(C(:)'*P0(:)) );
cvx_end

% recovery
[I,J] = find(P0 > 1e-7);

% ***********************
% ***********************


% ***** our solver ******
% ***********************
U = FFW(problem,options);

%% prony extraction
options_prony.factorized = 1;
[x,a] = mvprony(U,nn,options_prony);

% display
clf, hold on;
% ground truth
scatter(x1(I),x2(J),'o');
% ffw + prony
scatter(x(:,1),x(:,2),50,'x','linewidth',2);
scatter(x1,zeros(s,1),100,'.');
scatter(zeros(s,1),x2,100,'.');
