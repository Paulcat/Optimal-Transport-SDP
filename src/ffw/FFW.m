function [U,info] = FFW(problem,options)
%FFW Fourier-based Frank Wolfe algorithm
%   FFW(PROBLEM,OPTIONS) solves the user-specified SDP problem under Toeplitz constraint



% load problem parameters
m  = problem.vardim;
M  = prod(m);
d  = numel(m);
%
f   	= problem.fobj;
f0 	= problem.f0;
g   	= problem.grad;
%gt  	= problem.grad_pre;
lco 	= problem.ls;
type	= problem.name;
la 	= getoptions(problem,'hyper',0);
%scale = min(problem.hparams); % min la,rho to "normalize" criterion...
%
options_lmo  = set_lmo_options (options);
options_bfgs = set_bfgs_options(options);

pflag = strcmp(type,'Invariant');

% set options
U0    	= getoptions(options,'init',zeros(M,1)); %TODO: add error when initialization do not satisfy constraint...
Om    	= getoptions(options,'Om',ones(M,1));
maxit 	= getoptions(options,'maxiter',20);
tol   	= getoptions(options,'tol',1e-5);
verbose 	= getoptions(options,'display','on');
rho 		= getoptions(options,'rho',1e-3);


% Toeplitz penalization helper
[Tpen,Tpen_g,Tproj,~,Dnumel] = ffw_Tpen(m);


% *** Initialization ***
niter = 0;
crit  = -Inf;
U     = U0;
T 		= Tproj(U);
fval  = f(T) + 1/rho * Tpen(U,T);
v0    = ones(M,1)/sqrt(M);
%

% DISPLAY
if strcmp(verbose,'on')
	fprintf('\n\n')
	fprintf('-------------------------------- FFW Algorithm --------------------------------\n');
	fprintf('%-25s %-25s %-25s\n', 'GENERAL OPTIONS', 'LMO OPTIONS', 'BFGS OPTIONS');
	fprintf('%10s = %-11i  %10s = %-11i  %10s = %-11i\n', ...
		'maxiter', maxit, 'maxiter', options_lmo.maxiter, 'maxiter', options_bfgs.MaxIter);
	fprintf('%10s = %-12.1e %10s = %-12.1e %10s = %-12.1e\n', ...	
		'tol', tol, 'tol', options_lmo.tol, 'tol', options_bfgs.optTol);
	fprintf('%10s = %-10.4e\n', 'lambda', la);
	fprintf('%10s = %-10.4e\n', 'rho', rho); 
	fprintf('-------------------------------------------------------------------------------\n');
	fprintf('%-3s %-18s %-14s %-5s %-8s %-12s %-5s %-8s', ...
		'IT','OBJ','GAP','PI','(TIME)','LS:NEW/OLD','BFGS','(TIME)');
	fprintf('\n');
	fprintf('-------------------------------------------------------------------------------\n');
	%fprintf('%-3i %-+4.4e\n',niter,f(T) + 1/rho*Tpen(U,T))
end




E = []; % objective values
while crit < -tol && niter < maxit

	if strcmp(verbose,'debug')
		checkgradient(@(U)f(Tproj(U)), @(U)2*g(Tproj(U),U),U);
		checkgradient(@(U)1/rho*Tpen(U,Tproj(U)), @(U)2/rho*Tpen_g(U,Tproj(U),U),U);
		checkgradient(@(U)f(Tproj(U)) + 1/rho*Tpen(U,Tproj(U)), @(U)2*(g(Tproj(U),U) + 1/rho*Tpen_g(U,Tproj(U),U)), U);
	end
	E = [E; fval];

	% *** Linear Minimization Oracle ***
	g_lmo  = @(h) g(T,h) + 1/rho * Tpen_g(U,T,h);
	g_lmo1 = @(h) Om .* g_lmo( Om.* h);
	tic;
	[eVecm,nPI] = ffw_lmo(g_lmo1,v0,options_lmo);
	time_lmo = toc;
	eVecm = Om .* eVecm;

	% stopping criterion
	ege  = eVecm' * g_lmo(eVecm);
	crit = ege;
	%crit = la*rho/(la+rho) * crit; % scaling?


	% *** Frank-Wolfe update (with linesearch) ***
	t = Tproj(eVecm);
	%co = lco(U,T,eVecm,t);
	[mu,nu,stop] = ffw_ls(lco,U,T,eVecm,t,f0,rho,Dnumel,pflag);
	if mu==0
		U = sqrt(nu)*eVecm;
	elseif nu==0
		U = sqrt(mu)*U;
	else
		U = [sqrt(mu)*U, sqrt(nu)*eVecm];
	end


	% *** BFGS step ***
	%profile on;
	if options_bfgs.on
		tic;
		[U,nBFGS] = ffw_bfgs(f,g,Tpen,Tpen_g,Tproj,U,rho,pflag,options_bfgs);
		time_bfgs = toc;
	end
	%profile viewer;

	% update
	niter = niter+1;
	T = Tproj(U);
	fval = f(T) + 1/rho * Tpen(U,T);

	% display
	if strcmp(verbose,'on')
		fprintf('%-3i %-+17.10e  %-+13.6e  %-5i (%5.1f)  %-11.4e  %-5i (%5.1f)\n', ...i
			niter,fval,crit,nPI,time_lmo,mu/nu,nBFGS,time_bfgs);
	end
end

if strcmp(verbose,'on')
	fprintf('-------------------------------------------------------------------------------\n');
	if crit >= -tol
		fprintf('iterations stopped: tolerance reached\n');
	else
		fprintf('maximum iterations reached\n');
	end
	fprintf('\n\n');
end

info.E 	 = E;
info.crit = crit; % final value of criterion

end



function opt_bfgs = set_bfgs_options(options)
opt_bfgs.on              = getoptions(options, 'bfgsOn', 1);
opt_bfgs.display         = 'off';                                    % off | final | iter | full | excessive
opt_bfgs.optTol          = 1e-16;                                    % first-order optimality (exitflag 1)
opt_bfgs.progTol         = getoptions(options, 'bfgsProgTol', 1e-8); % parameters change (exitflag 2)
opt_bfgs.MaxFunEvals     = 100000;
opt_bfgs.MaxIter         = getoptions(options, 'bfgsMaxIter', 500);  % (exitflag 0)
opt_bfgs.Method          = 'lbfgs';
opt_bfgs.DerivativeCheck = 'off';
opt_bfgs.Corr            = 15;
opt_bfgs.Damped          = 0;
opt_bfgs.numDiff         = 0; % use-provided gradient
%
opt_bfgs.reg             = getoptions(options, 'bfgsReg', Inf);
end

function opt_lmo = set_lmo_options(options)
opt_lmo.tol     = getoptions(options, 'lmoTol', 1e-16);
opt_lmo.maxiter = getoptions(options, 'lmoMaxIter', 1000); 
end

