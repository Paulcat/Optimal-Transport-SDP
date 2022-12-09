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
			tau = getoptions(options,'reg',inf);
			%
			F   = @(X,TX) f(TX) + 1/2/tau * (norm(X,'fro')^2-1)^2 ... % trace penalization TODO: 1 or M?
				+ 1/rho * Tpen(X,TX);
			G   = @(X,TX) 2*g(TX,X) + 2/tau * (norm(X,'fro')^2-1)*X ...
				+ 2/rho * Tpen_g(X,TX,X);

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
			problem.M = spherecomplexfactory(M,r);

		elseif strcmp(cflag,'none')
			problem.M = euclideancomplexfactory(M,r);

		else
			error('This constraint cannot be handled: %s', cflag);
		end

		F = @(X,TX) f(TX) + 1/rho * Tpen(X,TX);
		G = @(X,TX) 2*g(TX,X) + 2/rho * Tpen_g(X,TX,X);
		%
		problem.cost  = @(X) F( X, Tprojn(m,X) );
		problem.egrad = @(X) G( X, Tprojn(m,X) );

		%checkgradient(problem);

		[U,fU,info] = rlbfgs(problem,U0,options);
		
		niter = [info.iter];     infos.niter = niter(end);
		time  = [info.time];     infos.time  = time(end);
		val   = [info.cost];     infos.val   = val(end);
		ng    = [info.gradnorm]; infos.ng    = ng(end);
end

end
