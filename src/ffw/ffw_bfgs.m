function [U,nit] = ffw_bfgs(f_user,g_user,Tpen,Tpen_g,Tproj,U0,rho,flag_constraint,options)
%FFW_BFGS BFGS step in FFW algorithm

deal2 = @(varargin) deal(varargin{1:nargout});

M = size(U0,1);
r = size(U0,2);

% from complex to real and 'vice-versa'
reshc = @(u) reshape( u(1:length(u)/2)    , [length(u)/2/r, r] ) + ...
        1i * reshape( u(length(u)/2+1:end), [length(u)/2/r, r] );
%
flatc = @(Z) [real(Z(:)); imag(Z(:))];

if flag_constraint
	tau = getoptions(options,'reg',inf);
	f = @(X,TX) f_user(TX) + 1/2/tau * (norm(U,'fro')^2 - M)^2 ... % regularization for trace constraint
		+ 1/rho * Tpen(X,TX);
	g = @(X,TX) 2*g_user(TX) + 2/tau * (norm(U,'fro')^2 - M)*U ...
		+ 2/rho * Tpen_g(X,TX);
else
	f = @(X,TX) f_user(TX) 	 + 1/rho * Tpen(X,TX);
	g = @(X,TX) 2*g_user(TX,X) + 2/rho * Tpen_g(X,TX,X);
end

fg  = @(X,TX) deal2(f(X,TX), flatc(g(X,TX)));
fg1 = @(X) fg( X, Tproj(X) );
fg2 = @(Z) fg1(reshc(Z)); %TODO: is it useful to decompose this much?

%checkgradient(fb,gb,flatc(U));
%drawnow;
%pause(1);

[U,val,exitflag,output] = lbfgs(fg2,flatc(U0),options);
U = reshc(U);

nit = output.iterations;
end

