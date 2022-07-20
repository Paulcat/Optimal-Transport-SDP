function c = genmom(n,meas)
% GENMOM compute trigonometric moments of Radon measures


d = length(n); % dimension

x = meas.x; % support
a = meas.a; % amplitude

K  = compute_grid(n,'spectral-sym');
Kr = reshape(K,prod(2*n+1),1,d);
xr = reshape(x,1,size(x,1),d);

Fx = exp(-2i*pi * sum(Kr.*xr,3)); % Fourier matrix on the support
c 	= Fx * a(:); % moments
c  = ifftshift(reshape(c,2*n+1,1)); % 0th moment on top left corner

end
