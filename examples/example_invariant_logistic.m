% logistic map example for invariant measure problem

clear all

% ***** set up problem ****
% *************************
% generate logistic map (Lagrange discretization)
r = 2; %2 : single point, 4: fully chaotic, 3.6: intermediate
f = @(x) r*x.*(1-x);
Z = rand(1e5,1); % initial states
% Lagrangian discretization
clf;
for it = 1:50
	hist(Z,linspace(0,1,100)); axis tight;
	drawnow;
	Z = f(Z);
end

L = 500;
t = (0:L-1)'/L;
% constant interpolation of push-forward map
[~,id_f]	= min(abs(f(t)-t'),[],2);
f_interp = sparse(id_f, 1:L, ones(L,1), L, L);
% piecewise linear interpolation
%I = floor(f(x)*L); r = f(x)*L - I;
%M = sparse(1+[I;min(I+1,L)], [1:L+1, 1:L+1]', [1-r;r], L+1, L+1);

%%

% eigenvalues of push-forward
% ***************************
[eVec,eVal] = eig(full(f_interp)); eVal = diag(eVal);
[eVal_Mo,I] = sort(abs(eVal-1),'ascend'); % eigenvalues modulus? 1 (ie |f_interp(x)| = |x|?) in first place
eVal = eVal(I); eVec = eVec(:,I);
%
J = find(eVal_Mo <= .05); % extract eigenvalue close to 1
plot(t,real(eVec(:,J)),'linewidth',3);

%%
% eigenvalues of Fourier transform of push-forward?
% *************************************************
n  = 30;
Ft = exp(-2i*pi*(-n:n)'*t(:)');
A  = 1/L * Ft * f_interp * Ft';

% % simple check
% B = fft2( full(f_interp .* ( exp(2i*pi*n*t(:)) .* exp(2i*pi*n*t(:)') )) );
% B = 1/L * B(1:2*n+1,1:2*n+1);
% norm(A-B)

%
[eVec_F,eVal_F] = eig(A); eVal_F = diag(eVal_F);
[eVal_F_Mo,I_F] = sort(abs(eVal_F-1),'ascend');
eVal_F = eVal_F(I_F); eVec_F = eVec_F(:,I_F);
%
p = real(Ft' * eVec_F(:,1)); %?

% display
clf, plot(t,p,'linewidth',3);

%%
% hyper parameters
la  = 1e-3; % scaling???
rho = 1e-3;


% variable size and functionals
m = n+1;
B = ifftshift(A - eye(2*n+1));
[f,f0] = inv_fobj(m,B,la);
g = inv_fgrad(m,B,f0,la);

% load problem
problem = struct;
%
problem.name   = 'invariant';
problem.vardim = m;
problem.fobj   = f;
problem.f0     = f0;
problem.grad   = g;
problem.cflag  = 'trace'; % trace constraint from invariant pb formulation
problem.hyper  = la;
problem.ls     = inv_lscoeffs(m,B,f0,la);
% *************************
% *************************



% **** load solver options ***
% ****************************
options = struct;
%
options.tol 			= 1e-5;
options.init 			= ones(m,1)/sqrt(m);
%options.init			= [1;zeros(m-1,1)];
options.maxiter 		= 40;
options.bfgsToolbox  = 'manopt';
options.bfgsProgTol 	= 1e-16;
options.bfgsMaxIter 	= 500;
options.bfgsReg 		= 0; % regularization of constraint in bfgs. Scaling?
options.lmoTol 		= 1e-10;
options.lmoMaxIter 	= 1e3;
% ***********************
% ***********************




% **** cvx ground truth ***
% *************************
cvx_begin sdp
cvx_precision high
variable T0(m,m) hermitian
variable c0(2*n+1) complex
%
T0 >= 0
c0(m) == 1
for j=-n:n
	diag(T0,j) == c0(j+m);
end
minimize( 1/2 * pow_pos(norm(c0,'fro'),2) + 1/2/la*pow_pos(norm(A*c0-c0,'fro'),2) );
cvx_end

% Prony extraction
options_prony.factorized = 0;
[x0,a0] = mvprony(T0,n,options_prony);

% ************************
% ************************


% ******* our solver *****
% ************************
U = FFW(problem,options);

% Prony extraction
options_prony.factorized = 1;
[x,a] = mvprony(U,n,options_prony);

%% display
clf, hold on;
stem(x0,a0,'filled','linewidth',2);
stem(1-x,a,'filled','linewidth',2); % why 1-x????
warning('something fishy here, 1-x instead of x is weird');

