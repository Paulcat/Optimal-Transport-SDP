% Test multivariate Prony 2d

% generate a curve
d = 2;
s = 5000;

ms = 10*s;

name = 'algebraic-circle';

switch name
	case 'circle'
		t  = linspace(0,1,s);
		r  = .3;
		z  = .5+.5i + r*exp(2i*pi*t);
		a0 = ones(s,1)/s;

	case 'treffle'
		t  = linspace(0,1,s);
		z  = .15*(1.3+cos(12*pi*t)) .* exp(2i*pi*t);
		z  = z/max(abs(z)) * .35;
		z  = z + .5+.5i;
		a0 = ones(s,1)/s; % amplitudes

	case 'algebraic-circle'
		w  = @(t) acos(-1 + 2*5/8 ./ (cos(2*pi*t) + 1))/2/pi;
		q  = zeros(s,1);
		s4 = s/4;

		filename = ['points_curve_M',int2str(s),'_quarter.txt'];
		for i=1:4
			q( ((i-1)*s4+1):(i*s4) ) = load(filename)/2/pi;
		end
		
		x  = repsym_curve(q,w);
		z  = x(:,1)' + 1i*x(:,2)';
		z  = .5 + .5i + z;
		a0 = ones(s,1)/s;

		mo = load('moments_curve_n1000.txt');
		mvec = @(kv) mo(sub2ind([1001 1001],abs(kv(:,1))+1,abs(kv(:,2))+1));
end

x0 = [real(z)', imag(z)'];


% display
clf;
scatter(x0(:,1),x0(:,2),ms*a0,'.');
xlim([0,1]), ylim([0,1]);
axis equal


%% 

nvec = 10*ones(1,d);
N = prod(nvec+1);
fX = compute_grid(nvec,'spectral');
fX = reshape(fX,N,d);

F0 = exp(-2i*pi * (fX(:,1)*x0(:,1)' + fX(:,2)*x0(:,2)'));
U0 = sqrt(a0)' .* F0; % moment matrix factorization

% short test
if strcmp(name,'algebraic-circle')
	fY = compute_grid(nvec,'spectral-sym');
	fY = reshape(fY,[prod(2*nvec+1),d]);
	c  = mvec(fY);
	c  = reshape(c,[2*nvec+1,1]);
	c  = ifftshift(c);

	options_P.mom = c;
end

% adding some noise?


% prony recovery
options_P.tol = 1e-20;
options_P.jdiag = 'cardoso';
options_P.factorized = 1;
options_P.verbose = false;

[x,a,info] = mvprony(U0,nvec,options_P);
a = a/sum(a,1); % hack?

% display
clf, hold on;
scatter(x0(:,1),x0(:,2),ms*a0,'.');
scatter(x(:,1),x(:,2),ms*a,'.');

%% Wasserstein distance with respect to n
nmin = 1;
nmax = 100;
step = 2;

options_W.maxit = 100;
options_W.wasserstein = 2;
options_W.norm = 2;
options_W.progtol = 1e-12;
options_W.verbose = 'off';

W = zeros((nmax-nmin+1)/2,1);

k = 0;
for i=nmin:step:nmax
	i
	k = k+1;
	nvec = i*ones(1,d);
	N = prod(nvec+1);
	fX = compute_grid(nvec,'spectral');
	fX = reshape(fX,N,d);

	F0 = exp(-2i*pi * (fX(:,1)*x0(:,1)' + fX(:,2)*x0(:,2)'));
	U0 = sqrt(a0)' .* F0;

	[x,a] = mvprony(U0,nvec,options_P);
	a = a/sum(a,1); % HACK?

	OT   = W_Semidiscrete(x0,a0,x,a,options_W);
	W(k) = OT(end);
end



