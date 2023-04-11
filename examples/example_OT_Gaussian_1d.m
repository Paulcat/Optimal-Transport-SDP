% example of FFW algorithm for OT between gaussian measures

clear all


% ***** set up problem ******
% ***************************
% create marginals
s = 500;
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
[n1,n2] = deal(15,15);
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
rho = 1e-4*Cr; % toeplitz penalization


% transport cost (in Fourier domain)
cost = zeros(2*nn+1);
cost(n1:n1+2,n2:n2+2) = [0 0 -1/4; 0 1/2 0; -1/4 0 0];
% shift
cost = ifftshift(cost);



% variable's size and functionals
mm     = nn+1;
[f,f0] = ot_fobj(mm,cost,c1,c2,la);
g = ot_fgrad(mm,cost,c1,c2,f0,la);
Tpen = ffw_Tpen(mm);

% load problem
problem = struct;
%
problem.name 	= 'optimal-transport';
problem.vardim = mm;
problem.fobj   = f;
problem.f0     = f0;
problem.grad   = g;
%problem.grad_pre 	= gU;
problem.cflag  = 'none';
problem.hyper 	= la;
problem.ls 	   = ot_lscoeffs(mm,cost,c1,c2,f0,la);
% ***************
% ***************

% **** load solver options *****
% ******************************
options = struct;
%
options.tol 			= 1e-5;
options.maxiter 		= 50;
options.bfgsToolbox  = 'minfunc';
options.bfgsProgTol 	= 1e-16;
options.bfgsMaxIter 	= 500;
options.lmoTol 		= 1e-10;
options.lmoMaxIter 	= 1e3;
options.rho				= rho;
options.display 		= 'on';
% ***********************
% ***********************



% **** cvx ground truth ****
% **************************
% periodic cost
C = sin(pi*(x1(:)-x2(:)')).^2;

% cvx solver
cvx_begin
variable P0(s,s) nonnegative
%
sum(P0,1) == a2.';
sum(P0,2) == a1;
%
minimize( trace(C(:)'*P0(:)) );
cvx_end
val_cvx = trace(C(:)'*P0(:));

% recovery
[I,J] = find(P0 > 1e-7);
x0 = [x1(I),x2(J)];
a0 = P0(sub2ind([s,s],I,J));
%
[fY,fX] = meshgrid(0:n2,0:n1);
F0 = exp(-2i*pi * (fX(:)*x0(:,1)' + fY(:)*x0(:,2)'));
U0 = F0 .* sqrt(a0');
T0 = Tprojn(mm,U0);
%
c01 = ifftshift(exp(-2i*pi*(-n1:n1)'*x0(:,1)') * a0(:));
c02 = ifftshift(exp(-2i*pi*(-n2:n2)'*x0(:,2)') * a0(:));

% ***********************
% ***********************


% ***** our solver ******
% ***********************
[U,info] = FFW(problem,options);

% prony extraction
options_prony.factorized = 1;
options_prony.jdiag = 'cardoso';
[x,a] = mvprony(U,nn,options_prony);



% objective values
val_ffw = info.E(end);
val_gt_in_penalized = f(T0) + 1/rho*Tpen(U0,T0);
%
CC = sin(pi*(x(:,1)-x(:,2))).^2;
val_ffw_in_constrained = sum(CC(:).*a(:));

% marginals
m1 = ifftshift(exp(-2i*pi*(-n1:n1)'*x(:,1)') * a(:));
m2 = ifftshift(exp(-2i*pi*(-n1:n1)'*x(:,2)') * a(:));
err_m1 = norm(m1-c1,'fro')/norm(c1,'fro');
err_m2 = norm(m2-c2,'fro')/norm(c2,'fro');

% toeplitz
err_T = sqrt(2*Tpen(U,Tprojn(mm,U))) / norm(U'*U,'fro');


fprintf('INFOS\n');
fprintf('%35s : %-10.7f\n', 'FFW final value', val_ffw);
fprintf('%35s : %-10.7f\n', 'GT (CVX) solution in penalized obj', val_gt_in_penalized);
fprintf('%35s : %-10.7f\n', 'FFW solution in constrained obj', val_ffw_in_constrained);
fprintf('%35s : %-10.7f\n', 'CVX final value', val_cvx);
fprintf('%35s : %-10.7f\n', '1st marginal constraint violation', err_m1);
fprintf('%35s : %-10.7f\n', '2nd marginal constraint violation', err_m2);
fprintf('%35s : %-10.7f\n', 'Toeplitz constraint violation', err_T);



% display
clf, hold on;
% ground truth
scatter(x1(I),x2(J),'o');
% ffw + prony
scatter(x(:,1),x(:,2),50,'x','linewidth',2);
scatter(x1,zeros(s,1),100,'.');
scatter(zeros(s,1),x2,100,'.');
