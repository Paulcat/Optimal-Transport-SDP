function B = fwd_op(n,type,params)
% FWD_OP Forward operator for invariant measure problems
%    Inputs:
%    - f: 
%    - spectral cutoff
%    - type: logistic | circle
%    - params: parameters, specific to each type

L = 500; % discretization in spatial domain

switch type
	case 'logistic'
		t = (0:L-1)'/L;

		r = getoptions(params,'potential',3.6);

		f = @(x) r * x.*(1-x);

		% constant interpolation of push-forward map
		[~,id_f] = min(abs(f(t)-t'),[],2);
		f_interp = sparse(id_f, 1:L, ones(L,1), L, L);

		% % piecewise linear interpolation
		% I = floor(f(t)*L); r = f(t)*L - I;
		% f_interp = sparse(1+[I;min(I+1,L)], [1:L+1, 1:L+1]', [1-r;r], L+1, L+1);

		% in Fourier domain
		% from spatial domain
		Ft = exp(-2i*pi*(-n:n)'*t(:)');
		A  = 1/L * Ft * f_interp *Ft';

		% analytic expression
		M1 = (-n:n) ./ (-n:n)';
		M1(isnan(M1)) = 1;
		M1(abs(M1)==Inf) = 0;
		M  = (-n:n)' .* (1 + 1/r * M1).^2;
		B  = 1/4/sqrt(r) * (1-1i) * erfi((1+i)*sqrt(pi*r)) * exp(-1i/2*pi*r * M);

end
