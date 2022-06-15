% Test for OT ffw functionals...

s = 500;
t = linspace(0,1,s);

[n1,n2] = deal(15,15);
nn 	  = [n1,n2];


% first measure
mu1 = .001 + .3 * exp(-(t-.3).^2 ./ (.02)^2 ) + .5 * exp(-(t-.6).^2 ./ (.08)^2 );
mu1 = mu1/norm(mu1,1);
%
%M1 = exp(-2i*pi*(0:n)'*t(:)') .* sqrt(mu1);
%M1 = M1*M1'; % moment matrix
c1 = exp(-2i*pi*(-n:n)'*t(:)') * mu1(:); % moment vector
c1 = ifftshift(c1);


% second measure
mu2 = .001 + .4 * exp(-(t-.5).^2 ./ (.04)^2 );
mu2 = mu2/norm(mu2,1);
%
%M2 = exp(-2i*pi*(0:n)'*t(:)') .* sqrt(mu2);
%M2 = M2*M2'; % moment matrix
c2 = exp(-2i*pi*(-n:n)'*t(:)') * mu2(:); % moment vector
c2 = ifftshift(c2);


% cost
h = zeros(2*nn+1);
h(n1:n1+2, n2:n2+2) = [0 0 -1/4; 0 1/2 0; -1/4 0 0];
h = ifftshift(h);


% parameters
[la,rho] = deal(rand,rand);


f = ot1_fobj (nn+1,h,c1,c2,la,rho); % TODO: why do I force n+1 here?
g = ot1_fgrad(nn+1,h,c1,c2,la,rho);
g = @(u) 2*g(u,u);
U = rand(prod(nn+1),6) + 1i*rand(prod(nn+1),6);
checkgradient(f,g,U);

%clf, hold on;
%area(t,mu1);
%area(t,mu2);



