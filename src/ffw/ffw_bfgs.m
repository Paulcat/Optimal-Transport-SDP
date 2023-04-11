function [U,infos] = ffw_bfgs(f,g,Tpen,Tpen_g,m,U0,rho,cflag,options)
%FFW_BFGS BFGS step in FFW algorithm

deal2 = @(varargin) deal(varargin{1:nargout});

M = size(U0,1);
r = size(U0,2);

solver = getoptions(options,'toolbox','minfunc');


switch solver
	case 'minfunc'
		% need to reshape complex matrices into real ones, and reciprocally
		flatc = @(Z) [real(Z(:)); imag(Z(:))];
		reshc = @(u) reshape( u(1:length(u)/2), [length(u)/2/r, r] ) + ...
			1i * reshape( u(length(u)/2+1:end), [length(u)/2/r, r] );

		% penalize potential aditional constraints (Moreau envelope)
		if strcmp(cflag,'trace') 
			tau = getoptions(options,'trace_reg',0);
			%
			if tau==0
				error('the bfgs of minfunc cannot handle constraints');
			else
				F = @(X,TX) f(TX) + 1/rho * Tpen(X,TX) + ...
					1/2/tau * (norm(X,'fro')^2-M)^2; % trace penalization
				G = @(X,TX) 2*g(TX,X) + 2/rho * Tpen_g(X,TX,X) + ...
					2/tau * (norm(X,'fro')^2-M)*X;
			end

		elseif strcmp(cflag,'none')
			F = @(X,TX) f(TX) + 1/rho * Tpen(X,TX);
			G = @(X,TX) 2*g(TX,X) + 2/rho * Tpen_g(X,TX,X);

		else
			error('This constraint cannot be handled: %s', cflag);
		end

		FG  = @(X,TX) deal2(F(X,TX), flatc(G(X,TX)));
		FG1 = @(X) FG( X, Tprojn(m,X) );
		FG2 = @(Z) FG1(reshc(Z)); %TODO: is it slower to decompose this much?

		% TODO: check gradient
		%checkgradient(fb,gb,flatc(U0));
		%drawnow

		tic;
		[U,val,exitflag,output] = lbfgs(FG2,flatc(U0),options);
		time = toc;

		U = reshc(U);

		infos.niter = output.iterations;
		infos.ng    = output.firstorderopt;
		vals        = output.trace.fval;
		infos.val   = vals(end);
		infos.time  = time;



	case 'manopt'
		if strcmp(cflag,'trace')
			tau = getoptions(options,'trace_reg',0);

			if tau==0
				problem.M = spherecomplexfactory(M,r); % trace of matrix = 1
				%
				F = @(X,TX) f(TX) + 1/rho * Tpen(X,TX);
				G = @(X,TX) 2*g(TX,X) + 2/rho * Tpen_g(X,TX,X); % euclidean gradient wrt U

			else
				warning('manopt can handle the trace constraint without relaxation');

				problem.M = euclideancomplexfactory(M,r);
				%
				F = @(X,TX) f(TX) + 1/rho * Tpen(X,TX) + ...
					1/2/tau * (norm(X,'fro')^2-1)^2;
				G = @(X,TX) 2*g(TX,X) + 2/rho * Tpen_g(X,TX,X) + ...
					2/tau * (norm(X,'fro')^2-1)*X;

			end

		elseif strcmp(cflag,'none')
			warning('it is not necessary to use manopt in the unconstrained case');

			problem.M = euclideancomplexfactory(M,r);
			%
			F = @(X,TX) f(TX) + 1/rho * Tpen(X,TX);
			G = @(X,TX) 2*g(TX,X) + 2/rho * Tpen_g(X,TX,X);

		else
			error('This constraint cannot be handled: %s', cflag);
		end

		problem.cost  = @(X) F( X, Tprojn(m,X) );
		problem.egrad = @(X) G( X, Tprojn(m,X) );

		%checkgradient(problem);

		[U,fU,info] = rlbfgs(problem,U0,options);

		if strcmp(cflag,'trace')
			U  = sqrt(M) * U; % rescaling so that trace(U*U') = M
			fU = F(U,Tprojn(m,U)); % fU = M^2 * fU would suffice for the invariant functional
		end
		
		infos.niter = info(end).iter;
		infos.time  = info(end).time;
		%infos.val  = info(end).cost;
		infos.val   = fU;
		infos.ng    = info(end).gradnorm;
end

end
