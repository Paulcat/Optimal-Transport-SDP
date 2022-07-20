% Benchmark algorithm for OT



% **** General set up **** %
d = 1; % dimension of marginals
[s1,s2] = deal(4,3); % sparsity of marginals

% 1st marginal
mu1.type			= 'discrete';
mu1.sparsity	= s1;
%
mu1 = genmeas(d,mu1);

% 2nd marginal
mu2.type			= 'discrete';
mu2.sparsity	= s2;
%
mu2 = genmeas(d,mu2);

% cost function: sin^2(x-y)
costF = [0 0 -1/4; 0 1/2 0; -1/4 0 0]; % in Fourier
cdegr	= 1; % only (0,0), (1,-1) and (-1,1) are non-zero
% ************************ % 




% **** BENCHMARK: HYPERPARAMETERS ****
scales_l = -8:2:8; L = length(scales_l);
scales_r = -8:2:8; R = length(scales_r);


% data = trigonometric moments
[n1,n2]	= deal(10,10);
if min(n1,n2) <= cdegr
	error('n should be at least greater than the degree of the cost function');
end
%
nn			= [n1,n2]; mm = nn+1;
c1			= genmom(n1,mu1);
c2			= genmom(n2,mu2);
costF		= padarray(costF,nn-1);



% scaling of hyperparameters
Cl = (norm(eval1(c1),'inf') + norm(eval2(c2),'inf'))/2;
Cr = prod(nn) * (norm(c1,'inf') + norm(c2,'inf'))/2;



% cvx "ground truth"
cost = sin(mu1.x(:) - mu2.x(:)').^2;

cvx_begin quiet
variable P0(s1,s2) nonnegative
% constraints
sum(P0,1) == a2.';
sum(P0,2) == a1;
%
minimize( trace(cost(:)'*P0(:)) );
cvx_end

[I,J] = find(P0 > 1e-7); % recovery
x0		= [mu1.x(I), mu2.x(J)];



% storage for benchmarks
fval	= zeros(L,R);
gnorm = zeros(L,R);
scrit = zeros(L,R);
errx	= zeros(L,R);



% let's go
problem = struct;
problem.name	= 'OT';
problem.vardim = mm;
%parpool(4);
parfor i=1:L
	i
	la = Cl*10.^scales_l(i);
	
	for j=1:R
		rho = Cr*10.^scales_r(j);
		
		% functionals
		[f,f0] = ot1_fobj		(mm,costF,c1,c2,	 la,rho);
		[g,gU] = ot1_fgrad	(mm,costF,c1,c2,f0,la,rho);
		ls		 = ot1_lscoeffs(mm,costF,c1,c2,f0,la,rho);
		
		% set problem
		problem = struct;
		problem.name		= 'OT';
		problem.vardim		= mm;
		problem.fobj		= f;
		problem.grad		= g;
		problem.grad_pre	= gU;
		problem.ls			= ls;
		
		% set solver options
		options = struct;
		options.display		= 'off';
		options.tol				= 1e-5;
		options.maxiter		= 30;
		options.bfgsProgTol	= 1e-16;
		options.bfgsMaxIter	= 500;
		options.lmoTol			= 1e-10;
		options.lmoMaxIter	= 1e3;
		
		% **** solve ****
		[U,infos] = FFW(problem,options);
		%
		options_prony = struct;
		options_prony.factorized = 1;
		[x,a] = mvprony(U,nn,options_prony);
		
		% store
		normg = 1;
		for k=1:50
			v	= rand(prod(nn+1),1);
			v	= v/norm(v,'fro');
			normg = max(normg, norm(g(U,v),'fro') / norm(v,'fro');
		end
		%
		fval (i,j) = f(U);
		gnorm(i,j) = normg;
		errx (i,j) = Jaccard(x0,x,.1);
		scrit(i,j) = infos.crit;
	end
end
% **************************************** %



%% display
clf;

subplot(2,2,1);
[B,A] = meshgrid(Cl*10.^scales_l,Cr*10.^scales_r);
surf(A,B,log10(gnorm)), view(2), colorbar;
set(gca,'xscale','log','yscale','log');
set(gca,'XLimSpec','tight','YlimSpec','tight');
title('Norm of gradient (log-scale)');
xlabel('\rho');
ylabel('\lambda');

subplot(2,2,2);
surf(A,B,log10(abs(scrit))), view(2), colorbar;
set(gca,'xscale','log','yscale','log');
set(gca,'XLimSpec','tight','YlimSpec','tight');
title('Final criterion (absolute value, log-scale)');
xlabel('\rho');
ylabel('\lambda');

subplot(2,2,3);
surf(A,B,log10(abs(fval))), view(2), colorbar;
set(gca,'xscale','log','yscale','log');
set(gca,'XLimSpec','tight','YlimSpec','tight');
title('Final value (absolute value, log-scale)');
xlabel('\rho');
ylabel('\lambda');

subplot(2,2,4);
surf(A,B,errx), view(2), colorbar;
set(gca,'xscale','log','yscale','log');
set(gca,'XLimSpec','tight','YlimSpec','tight');
title('Jaccard index');
xlabel('\rho');
ylabel('\lambda');
