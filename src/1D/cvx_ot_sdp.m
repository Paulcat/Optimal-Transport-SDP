function T = cvx_ot_sdp(n,cost,u1,u2,la,solver)
% CVX_OT_SDP Solve Lasserre's relaxation of OT with cvx


N = prod(n+1);
X = compute_grid(n,'spectral');
X = reshape(X,N,2); % frequency grid

ids = marginals(n,'colex',1); % linear indices of moments in matrix without repetition
i1  = ifftshift(ids{1}); % linear indices of 1st marginal moments in matrix
i2  = ifftshift(ids{2}); % linear indices of 2nd marginal moments in matrix


switch solver
	case 'primal'
		cvx_precision best
		cvx_solver SDPT3
		cvx_begin sdp
		variable T(N,N) hermitian % moment matrix
		variable c(2*n+1) complex % moments
		%
		R >= 0
		T = zeros(N);
		for k1=-n(1):n(1)
			for k2=-n(2):n(2)
				Th = kron( toepl(k2), toepl(k1) );
				T  = T + c(k1+n(1)+1, k2+n(2)+1)*Th;
			end
		end
		R == T;
		%

		if la == 0
			R(i1) == u1; % 1st marginal (moment) constraint
			R(i2) == u2; % 2nd marginal ...
			minimize( real(trace(cost'*ifftshift(c))) );
		else
			minimize( real(trace(cost'*ifftshift(c))) + ...
				1/4/la*norm(R(i1)-u1,'fro')^2 + 1/4/la*norm(R(i2)-u2,'fro')^2 );
		end
		cvx_end
end

end
