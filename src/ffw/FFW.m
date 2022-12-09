function [U,info] = FFW(problem,options)
%FFW Fourier-based Frank Wolfe algorithm
%   FFW(PROBLEM,OPTIONS) solves the user-specified SDP problem under Toeplitz constraint



% load problem parameters
m = problem.vardim;
M = prod(m);
d = numel(m);
%
f     = problem.fobj;
f0    = problem.f0;
g     = problem.grad;
lco   = problem.ls;
type  = problem.name;
cflag = problem.cflag;
la    = getoptions(problem,'hyper',0);
%scale = min(problem.hparams); % min la,rho to "normalize" criterion...
%
options_lmo  = set_lmo_options (options);
options_bfgs = set_bfgs_options(options);

% set options
U0       = getoptions(options,'init',zeros(M,1)); %TODO: add error when initialization do not satisfy constraint...
Om       = getoptions(options,'Om',ones(M,1));
maxit    = getoptions(options,'maxiter',20);
tol      = getoptions(options,'tol',1e-5);
verbose 	= getoptions(options,'display','on');
rho      = getoptions(options,'rho',1e-3);


% Toeplitz penalization helper
[Tpen,Tpen_g] = ffw_Tpen(m);


% *** Initialization ***
niter = 0;
crit  = -Inf;
U     = U0;
T 	   = Tprojn(m,U);
fval  = f(T) + 1/rho * Tpen(U,T);
v0    = ones(M,1)/sqrt(M);
%

% DISPLAY
if strcmp(verbose,'on')
   fprintf('\n\n')
   fprintf('---------------------------------------------- FFW Algorithm -------------------------------------------------\n');
   fprintf('Fourier Frank-Wolfe solver for low-rank semidefinite programming, see %s\n', '<a href="https://epubs.siam.org/doi/10.1137/19M124071X"> Catala, Duval, Peyré [2019] </a>');
   fprintf('**\n');
   fprintf('%-10s : %-s\n', 'problem', type);
   fprintf('%-10s : %-i\n', 'dimension', d);
   fprintf('%-10s : ', 'order'); fprintf('%-g ', m); fprintf('(size = %i x %i)\n', M, M); 
   fprintf('--------------------------------------------------------------------------------------------------------------\n');
   fprintf('%-25s %-25s %-25s\n', '*GENERAL OPTIONS*', '*LMO OPTIONS*', '*BFGS OPTIONS*');
   fprintf('%10s = %-11i  %10s = %-11i  %10s = %-11i\n', ...
		'maxiter', maxit, 'maxiter', options_lmo.maxiter, 'maxiter', options_bfgs.MaxIter);
   fprintf('%10s = %-12.1e %10s = %-12.1e %10s = %-12.1e\n', ...	
		'tol', tol, 'tol', options_lmo.tol, 'tol', options_bfgs.optTol);
   fprintf('%10s = %-10.4e\n', 'lambda', la);
   fprintf('%10s = %-10.4e\n', 'rho', rho); 
   fprintf('--------------------------------------------------------------------------------------------------------------\n');
   fprintf('%-3s %-18s %-14s %-14s %-14s %-5s %-8s %-12s %-5s %-8s', ...
		'IT','OBJ','CRITERION','GAP CERTIF','GRAD SNORM','PI','TIME(s)','LS:OLD/NEW','BFGS','TIME(s)');
   fprintf('\n');
   fprintf('------------------------------------------------------------------------------------------------------------------\n');
	%fprintf('%-3i %-+4.4e\n',niter,f(T) + 1/rho*Tpen(U,T))
end




E_list    = []; % objective values
ng_list   = [];
gap_list  = [];
ttot_bfgs = 0;
ttot_lmo  = 0;
%
while crit < -tol && niter < maxit

	if strcmp(verbose,'debug')
		Tp = @(U) Tprojn(m,U);
		checkgradient(@(U)f(Tp(U)), @(U)2*g(Tp(U),U),U);
		checkgradient(@(U)1/rho*Tpen(U,Tp(U)), @(U)2/rho*Tpen_g(U,Tp(U),U),U);
		checkgradient(@(U)f(Tp(U)) + 1/rho*Tpen(U,Tp(U)), @(U)2*(g(Tp(U),U) + 1/rho*Tpen_g(U,Tp(U),U)), U);
	end
	E_list = [E_list; fval];

	% *** Linear Minimization Oracle ***
	g_lmo  = @(h) g(T,h) + 1/rho * Tpen_g(U,T,h);
	g_lmo1 = @(h) Om .* g_lmo( Om.* h);
	%
	[eVecm,iPI] = ffw_lmo(g_lmo1,v0,options_lmo);
	ttot_lmo    = ttot_lmo + iPI.time;
	
	% output of LMO
	eVecm    = Om .* eVecm; % eigenvector fr minimal eigenvalue

	% upper bound on dual gap TODO: how to scale it?
	ege      = eVecm' * g_lmo(eVecm);
	crit     = ege;
	%crit = la*rho/(la+rho) * crit; % scaling?
	gap      = real(trace(U'*g_lmo(U)) - ege);
	gap_list = [gap_list,gap];

	% check spectral norm of gradient
	ng      = iPI.eVmax; % maximal absolute eigenvalue of gradient
	ng_list = [ng_list, ng];

	if crit >= -tol
		break;
	end


	% *** Frank-Wolfe update (with linesearch) ***
	t = Tprojn(m,eVecm);
	%co = lco(U,T,eVecm,t);
	[mu,nu,stop] = ffw_ls(lco,m,U,T,eVecm,t,f0,rho,cflag);
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
		[U,iBFGS] = ffw_bfgs(f,g,Tpen,Tpen_g,m,U,rho,cflag,options_bfgs);
		%
		ttot_bfgs = ttot_bfgs + iBFGS.time;
	end
	%profile viewer;
	

	% update
	niter = niter+1;
	T     = Tprojn(m,U);
	fval  = f(T) + 1/rho * Tpen(U,T);

	% display
	if strcmp(verbose,'on')
		fprintf('%-3i %-+18.10e %-+14.6e %-+14.6e %-14.6e %-5i %-8.2f %-12.4e %-5i %-8.2f\n', ...
			niter,fval,crit,gap,ng,iPI.niter,iPI.time,mu/nu,iBFGS.niter,iBFGS.time);
	end
end

if strcmp(verbose,'on')
	fprintf('----------------------------------------------------------------------------------------------------------\n');
	if crit >= -tol
		fprintf('Iterations stopped: tolerance reached\n');
	else
		fprintf('Iterations stopped: maximum number reached\n');
	end
	fprintf('%20s : %-10.5f\n', 'objective', fval);
	fprintf('%20s : %-10.2f\n', 'time spent in LMO', ttot_lmo);
	fprintf('%20s : %-10.2f\n', 'time spent in BFGS', ttot_bfgs);
	fprintf('\n\n');
end

info.E     = E_list;
info.gap   = gap_list;
info.ngrad = ng_list;
info.crit  = crit; % final value of criterion

end



function opt_bfgs = set_bfgs_options(options)
opt_bfgs.toolbox         = getoptions(options,'bfgsToolbox','minfunc'); % minfunc | manopt

% options for minfunc
opt_bfgs.on              = getoptions(options, 'bfgsOn', 1);
opt_bfgs.display         = 'off'; % off | final | iter | full | excessive
opt_bfgs.optTol          = getoptions(options, 'bfgsOptTol', 1e-16); % first-order optimality (exitflag 1)
opt_bfgs.progTol         = getoptions(options, 'bfgsProgTol', 1e-8); % parameters change (exitflag 2)
opt_bfgs.MaxFunEvals     = 100000;
opt_bfgs.MaxIter         = getoptions(options, 'bfgsMaxIter', 500); % (exitflag 0)
opt_bfgs.Method          = 'lbfgs';
opt_bfgs.DerivativeCheck = 'off';
opt_bfgs.Corr            = 15;
opt_bfgs.Damped          = 0;
opt_bfgs.numDiff         = 0; % use-provided gradient
opt_bfgs.reg             = getoptions(options, 'bfgsReg', Inf); % regularization in case of additional constraint

% options for manopt
opt_bfgs.maxiter          = getoptions(options, 'bfgsMaxIter', 500);
opt_bfgs.verbosity        = 0;
opt_bfgs.tolgradnorm      = getoptions(options, 'bfgsOptTol', 1e-16);
end

function opt_lmo = set_lmo_options(options)
opt_lmo.tol      = getoptions(options, 'lmoTol', 1e-16);
opt_lmo.maxiter  = getoptions(options, 'lmoMaxIter', 1000); 
end

