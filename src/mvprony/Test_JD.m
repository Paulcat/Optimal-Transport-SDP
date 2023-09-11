% Wasserstein distance between discrete approximation and measure
% See [Cardoso, Catala, Duval, Peyré]

%if ~exist('results-ajd','dir')
%	mkdir('results-ajd');
%end

normalize1 = @(v) v/norm(v,1);

d = 2;

type = 'treffle';
s    = 5000;

switch type
	case 'treffle'
		t = linspace(0,1,s);
		z = .15*(1.3 + cos(6*pi*t)) .* exp(2i*pi*t);
		z = .15*(1.3 + cos(12*pi*t)) .* exp(2i*pi*t);
		z = z/max(abs(z)) * .35;
		z = z + .5 + .5i;

		S = [real(z(:)),imag(z(:))];
		A = ones(s,1)/s;
end


% load options for semidiscrete optimization
options_W.maxit = 100;
options_W.parallelize = 0;
options_W.multiscale = 0;
options_W.wasserstein = 2;
options_W.norm = 2;
options_W.progtol = 1e-12;


% load options for mvprony
options_P.tol = 1e-3;
options_P.jdiag = 'cardoso';
options_P.factorized = 1;



step = 2;
nmin = 1;
nmax = step*50 + nmin - 1;



W = zeros(1,(nmax-nmin+1)/step);

for n=1:(nmax-nmin+1)/2
	nvec = (step*n + (nmin-1))*ones(1,d);
	N = prod(nvec+1);

	% moment vector
	fX = compute_grid(nvec,'spectral');
	fX = reshape(fX,[N,d]);
	%
	U = sqrt(A)' .* exp(-2i*pi * (fX(:,1)*S(:,1)' + fX(:,2)*S(:,2)') );

	[x,a,info] = mvprony(U,nvec,options_P);

	OT = W_Semidiscrete(S,A,x,a,options_W);
	W(n) = OT(end);
end

%filename = ['results-ajd/W',int2str(options_W.wasserstein),'-',int2str(d),'D_',type,'_n',...
%	int2str(nmin),'-',int2str(nmax)];
%fullname = [filename, '.csv'];
%if isfile(fullname)
%	csvwrite([filename,'-bis.csv']);
%else
%	csvwrite(fullname,W);
%end

