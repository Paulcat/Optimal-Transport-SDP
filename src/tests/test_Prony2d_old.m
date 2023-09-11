% Test multivariate Prony

%clear all

% generate a curve
%name = 'sparse-rand';
%name = '3-sparse';
%name = 'sparse1';
%name = 'sparse2';
%name = 'circle';
%name = '2-ellipses';
%name = 'rand-curve';
name = 'mix';
%name = 'mix-1d';
%name = 'rand-area';
%name = 'treffle';
%name = 'line1';
%name = 'line2';
%name = 'lissajou1';
%name = 'lissajou2';
%name = 'helix';
%name = 'disc';
%name = 'disc-random';
%name = 'pseudo-disc';
%name = 'square-random';

ms = 50;

d = 2;

switch name
	case '3-sparse'
		n = 3;
		l = 1e-8;
		x = [.4; .4+l/2; .4+l];
		y = [.4; .4+sqrt(3)/2*l; .4];
		z = x+1i*y;
		%a = rand(n,1) + 1i * rand(n,1);
		a = 0.5*[1; 1; 1];
		ms = 300;
		
	case 'circle'
		n = 5000;
		r = .3;
		z = .5+.5i + r*exp(2i*pi*linspace(0,1,n));
		%offset = r/sqrt(2);
		%z = z + offset*(1-1i);
		%z = z + (.2 - .2i);
		a = ones(n,1)/n; % amplitudes
		%a = rand(n,1) + 1i * rand(n,1); % amplitudes
		
		
	case '2-ellipses'
		n = 5000;
		t = linspace(0,1,n/2);
		
		% ellipse 1
		a1 = .4; b1 = .4;
		c1 = .6; d1 = .6;
		e1 = .7;
		al1 = 1/2 * sqrt((c1-a1)^2 + (d1-b1)^2);
		be1 = al1*(1-e1^2);
		X1 = al1*cos(2*pi*t);
		Y1 = be1*sin(2*pi*t);
		w1 = atan2(d1-b1,c1-a1);
		ell1x = (a1+c1)/2 + X1*cos(w1) - Y1*sin(w1);
		ell1y = (b1+d1)/2 + X1*sin(w1) + Y1*cos(w1);
		
		% ellipse 2
		a2 = .7; b2 = .3;
		c2 = .3; d2 = .7;
		e2 = .9;
		al2 = 1/2 * sqrt((c2-a2)^2 + (d2-b2)^2);
		be2 = al2*(1-e2^2);
		X2 = al2*cos(2*pi*t);
		Y2 = be2*sin(2*pi*t);
		w2 = atan2(d2-b2,c2-a2);
		ell2x = (a2+c2)/2 + X2*cos(w2) - Y2*sin(w2);
		ell2y = (b2+d2)/2 + X2*sin(w2) + Y2*cos(w2);
		
		z = [ell1x,ell2x] + 1i*[ell1y,ell2y];
		
		a = ones(n,1)/n;
		a(n/2+1:end) = 2*a(n/2+1:end);
		
	case 'rand-curve'
		n = 1000;
		
		% Randomize amplitude and phase.
		H = 20;
		rho = rand(1,H) .* logspace(-0.5,-2.5,H);
		phi = rand(1,H) .* 2*pi;
		
		% Accumulate r(t) over t=[0,2*pi]
		t = linspace(0,2*pi,n);
		r = ones(size(t));
		for h=1:H
			r = r + rho(h)*sin(h*t+phi(h));
		end
		
		% Reconstruct z(t)
		x = .5 + .2*r .* cos(t);
		y = .5 + .2*r .* sin(t);
		z = x + 1i*y;
		
		a = ones(n,1)/n;
		
	case 'rand-area'
		% Randomize amplitude and phase.
		H = 20;
		rho = rand(1,H) .* logspace(-0.5,-2.5,H);
		phi = rand(1,H) .* 2*pi;
		x = [];
		y = [];
		for i=1:100
			
			% Accumulate r(t) over t=[0,2*pi]
			t = linspace(0,2*pi,5*i);
			r = ones(size(t));
			for h=1:H
				r = r + rho(h)*sin(h*t+phi(h));
			end
			
			% fill up area
			x = [x, .5 + i/100*.2*(r+.15*rand(size(r))).*cos(t+.15*rand(size(t)))];
			y = [y, .5 + i/100*.2*(r+.15*rand(size(r))).*sin(t+.15*rand(size(t)))];
			%x = [x, .5 + i/100*.2*(r).*cos(t)];
			%y = [y, .5 + i/100*.2*(r).*sin(t)];
		end
		z = x(:) + 1i*y(:);
		a = ones(size(z))/length(z);
		
	case 'mix'
		n1 = 500;
		t1 = linspace(0,2*pi,n1);
		x1 = .3 + .15*cos(t1);
		y1 = .6 + .15*sin(t1);
		a1 = 3*ones(1,n1)/n1;
		
% 		n2 = 500;
% 		t2 = linspace(0,2*pi,4*n2);
% 		[Y,X] = meshgrid(t2);
% 		I  = find((X(:)-.6).^2+(Y(:)-.3).^2 <= .1^2);
% 		x2 = X(I)';
% 		y2 = Y(I)';
% 		a2 = ones(1,n2)/n2;

		n2 = 500;
		x2 = linspace(.7,.8,n2);
		y2 = mod(2*pi*x2/3+.8,1);
		a2 = ones(1,n2)/n2;

		
		n3 = 6;
		x3 = rand(1,n3);
		y3 = rand(1,n3);
		a3 = rand(1,n3);
		a3 = a3/sum(a3);
		
		n = n1+n2+n3;
		x = [x1,x2,x3];
		y = [y1,y2,y3];
		a = [a1,a2,a3];
		
		%load('example-misc-2.mat');
		%a(1:1000) = [a1,a2];
		%z = x+1i*y;
		
		
% 		n1 = 500;
% 		% Randomize amplitude and phase.
% 		H1 = 20;
% 		rho1 = rand(1,H1) .* logspace(-0.5,-2.5,H1);
% 		phi1 = rand(1,H1) .* 2*pi;
% 		
% 		% Accumulate r(t) over t=[0,2*pi]
% 		t1 = linspace(0,2*pi,n1);
% 		r1 = ones(size(t1));
% 		for h=1:H1
% 			r1 = r1 + 1.5*rho1(h)*sin(h*t1+phi1(h));
% 		end
% 		
% 		a1 = ones(n1,1)/n1;
% 		
% 		n2 = 9000;
% 		x2 = -1+2*rand(1,n2*500);
% 		y2 = -1+2*rand(1,n2*500);
% 		I = find(sqrt((x2.*cos(5*pi/6)+y2.*sin(5*pi/6)).^2/.06^2+...
% 			(x2.*sin(5*pi/6)+y2.*cos(5*pi/6)).^2/.08^2)<=1);
% 		x2 = .6 + x2(I(1:n2));
% 		y2 = .3 + y2(I(1:n2));
% 		
% 		a2 = ones(n2,1)/n2;
% 		
% 		n3 = 6;
% 		a3 = rand(n3,1);
% 		a3 = a3/sum(a3)/3;
% 		n  = n1+n2+n3;
% 		
% 		% Reconstruct z(t)
% 		x = [.4+.1*r1.*cos(t1), x2, rand(1,n3)];
% 		y = [.7+.1*r1.*sin(t1), y2, rand(1,n3)];
% 		z = x + 1i*y;
% 		
% 		a = [a1;a2;a3];
% 		a = a/sum(a);
% 		
% 		clear all
% 		name = 'mix';
% 		
		n1 = 500;
		n2 = 9000;
		n3 = 6;
		n = n1+n2+n3;
		load('example-misc.mat');
		x = [x;.6]; y=[y;.84];
		n3 = n3+1;
		n = n+1;
		z = x+1i*y;
		a(end+1) = a(end);
		a = a/norm(a,1);
		
	case 'mix-1d'
		nn=2048;
		
		%n=19;
		
		%xx=0:1/nn:1;
		%xx(end)=[];
		
		%mu=zeros(nn,1);
		%mu(2*nn/8:5*nn/8-1)=8/9;
		
		%x0=1/8;
		%mu(nn*x0)=nn/3;
		
		%xx1=(6*nn/8:8*nn/8-1)/nn;
		%mu(6*nn/8:8*nn/8-1)=(1./sqrt(abs(xx1-7/8))-sqrt(8))/(3*sqrt(2));
		%mu(7/8*nn)=mu(7/8*nn-1);
		
		x1 = 1/8;
		a1 = 1/3;
		a1 = a1/sum(a1);
		
		x2 = linspace(1/4,5/8,500);
		a2 = 8/9*ones(1,500)/500;
		a2 = a2/sum(a2);
		
		x3 = linspace(3/4,1,500);
		a3 = sqrt(2)/3*(1./sqrt(abs(x3-7/8)) - sqrt(8));
		a3 = a3/sum(a3);
		
		%muh=fft(mu);
		%k=0:nn-1;
		%k=min(k',nn-k');
		%Fnh=1-k/(n+1);
		%Fnh(Fnh<0)=0;
		%ph=muh.*Fnh;
		%p=ifft(ph);
		
		%Dnh=ones(size(muh));
		%Dnh(k>n)=0;
		%ph0=muh.*Dnh;
		%p0=ifft(ph0);
		
		%c=3/4;
		%mu(nn*x0)=0;
		
		x = [x1,x2,x3];
		a = [a1,a2,a3];
		a = a/sum(a);
		
		z = x + 1i*0;
		
	case 'treffle'
		% 1D curve
		n = 5000;
		t = linspace(0,1,n);
		z = .15*(1.3+cos(6*pi*t)) .* exp(2i*pi*t);
		z = .15*(1.3+cos(12*pi*t)) .* exp(2i*pi*t);
		z = z/max(abs(z)) * .35;
		z = z + .5+.5i;
		a = ones(n,1)/n; % amplitudes
		
	case 'lissajou1'
		n = 5000;
		t = linspace(0,1,n);
		z = 3/4 * cos(pi*t).^3 .* sin(pi*t) + 3i/4*(cos(pi*t).^2);
		z = z + .5+.1i;
		a = ones(n,1)/n;
		
	case 'lissajou2'
		n = 5000;
		t = linspace(0,1,n);
		z = 1/4*cos(2*pi*t) + 1/6*cos(6*pi*t) + 1i/12 * (3*sin(2*pi*t)-2*sin(6*pi*t));
		z = z + .5+.5i;
		a = ones(n,1)/n;
		%a = rand(n,1) + 1i * rand(n,1); % amplitudes
		
	case 'disc'
		n = 300;
		t = linspace(0,1,n);
		[Y,X] = meshgrid(t,t); z = X(:)+1i*Y(:);
		I = find(abs(z-(.5+1i*.5))<.25);
		z = z(I);
		a = ones(n,1)/n; % amplitudes
		
	case 'disc-random'
		n = 20000;
		radius = 0.2;
		xc = 0.5; yc = 0.5;
		th = rand(n,1)*(2*pi);
		r = sqrt(rand(n,1))*radius;
		x = xc + r.*cos(th);
		y = yc + r.*sin(th);
		
		z = x + 1i*y;
		a = ones(n,1)/n;
		
	case 'sparse-rand'
		n = 10;
		z = 0.8*rand(n,1) + 0.8i * rand(n,1);
		z = z + .1+.1i;
		a = -1 + 2*rand(n,1);
		%a = rand(n,1) + 1i * rand(n,1); % amplitudes
		a = ones(n,1)/n;
		
		ms = 300;
		
	case 'sparse1'
		n = 2;
		%z = [0.3 0.6] + 1i*[0.3 0.6];
		z = [0.3 0.3] + 1i * [0.3 0.6];
		a = ones(n,1);
		
	case 'sparse2'
		n = 2;
		%z = [0.3 0.6] + 1i * [0.3 0.6];
		th = angle(-.5 * exp(2i*pi*s0) + exp(2i*pi*t0));
		t1 = mod(th/pi-t0,1);
		s1 = mod(th/pi-s0,1);
		z = [t0, t1] + 1i*[s0, s1];
		a = ones(n,1);
		
	case 'line1'
		m = 500;
		t = linspace(0,1,m);
		z = t + 1i*t;
		%z = rand(m,1) + 1i*rand(m,1);
		a = ones(m,1)/m;
		
	case 'line2'
		m = 500;
		x = rand(m-2,1); x = [0;sort(x);1];
		y = rand(m-2,1); y = [0;sort(y);1];
		z = x+1i*y;
		
		a = ones(m,1)/m;
		
	case 'helix'
		n = 1000;
		t = linspace(0,1,n);
		z = 1/3 * cos(pi*t).*sin(pi*t) - 1i/3*sin(2*pi*t).*cos(2*pi*t);
		z = z + .5+.5i;
		a = ones(n,1)/n;
		
	case 'pseudo-disc'
		n = 1000;
		th = linspace(0,100*pi,n);
		r = linspace(0,0.2,n);
		z = r.*exp(1i*th);
		z = z + .5+.5i;
		a = ones(n,1)/n;
		%a = rand(n,1) + 1i * rand(n,1); % amplitudes
		
	case 'square-random'
		sizeI = [1024,1024];
		spacing = 10;
		pts = poissonDisc(sizeI,spacing);
		z = pts(:,1) + 1i*pts(:,2);
		z = z/max(z)/1.5;
		z = z + .18+.5i;
		%a = ones(length(z),1)/length(z);
		n = length(z);
		a = rand(n,1) + 1i * rand(n,1); % amplitudes
		
end
x = real(z)'; y = imag(z)';

clf;
if strcmp(name,'mix')
	hold on;
	%scatter(x(1:n1),y(1:n1),50,[0.8500 0.3250 0.0980],'.');
	%fill(x(n1+1:n1+n2),y(n1+1:n1+n2),[0.8500 0.3250 0.0980],'linestyle','none');
	%scatter(x(n1+1:n1+n2),y(n1+1:n1+n2),50,[0.8500 0.3250 0.0980],'.');
	%scatter(x(n1+n2+1:n),y(n1+n2+1:n),80,[0.8500 0.3250 0.0980],'filled');
	n3 = 6;
	scatter(x(1:end-n3),y(1:end-n3),50,[0.8500 0.3250 0.0980],'.');
	scatter(x(end-n3+1:end),y(end-n3+1:end),2e4*a(end-n3+1:end),[0.8500 0.3250 0.0980],'.');
	xlim([0,1]), ylim([0,1]);
	xticks([]), yticks([]);
	box on
	
	% plot polynomial
	fc = 13;
	n = [fc fc];
	N = 2*n+1;
	% frequency
	[fY,fX] = meshgrid(-fc:fc,-fc:fc);
	F = @(x,y)exp( -2i*pi*(fX(:)*x(:)' + fY(:)*y(:)') );
	% evaluate a polynomial associated to moment matrix on a grid
	T = 256;
	%t = linspace(0,1,m);
	t = (0:T-1)/T;
	[Y,X] = meshgrid(t,t); %X = X(:); Y = Y(:);
	%evalP = @(Mi)real( sum( (Mi*F(X,Y)) .* conj(F(X,Y)) , 1) );
	%evalP = @(Mi)reshape(evalP(Mi),[T,T]);
	evalP = @(y) reshape(real(F(X,Y)'*y(:)),[256,256]);
% 	
 	%U = sqrt(a).*F(x,y)';
 	%M = U'*U;
	surf(X,Y,evalP(F(x,y)*a(:)),'linestyle','none');
	%surf(X,Y,evalP(M),'linestyle','none');
	colormap(turbo);
	view(2);
	

elseif strcmp(name,'mix-1d')
	hold on
	%col = [0.8500 0.3250 0.0980];
	col = [.7 .7 .7];
	h = stem(x1,.005+a1,'.','color',col,'linewidth',2,'markersize',30);
	plot([1/4,x2,5/8],.005+70*[0,a2,0],'color',col,'linewidth',3);
	plot(x3,.005+70*a3,'color',col,'linewidth',3);
	plot([0,1/4],.005+[0 0],'color',col,'linewidth',2);
	plot([5/8,3/4],.005+[0 0],'color',col,'linewidth',2);
	xlim([0,1]);
	%axis tight
	h.ShowBaseLine = 'off';
	set(gca,'FontSize',10);
	xticks([]);
	yticks([]);
	box on;
 	%set(gca,'ylimspec','tight');
	xlim([0,1]), ylim([-.4476,3.113]);
	
	% plot polynomial
	clf, hold on;
	h = stem(x1,a1,'.','color',[.7 .7 .7],'linewidth',2,'markersize',30);
	h.ShowBaseLine = 'off';
	plot([1/4,x2,5/8],70*[0,a2,0],'color',[.7 .7 .7],'linewidth',3);
	plot(x3,70*a3,'color',[.7 .7 .7],'linewidth',3);
	plot([0,1/4],[0 0],'color',[.7 .7 .7],'linewidth',2);
	plot([5/8,3/4],[0 0],'color',[.7 .7 .7],'linewidth',2);
	fc = 13;
	N = 2*fc+1;
	% frequency
	F  = @(x) exp(-2i*pi*(-fc:fc)'*x(:)');
	%F  = @(x) exp(-2i*pi*(0:fc)'*x(:)');
% 	
	T = 256;
	t = (0:T-1)'/T;
	%evalP = @(Mi)1/T*real( sum( (Mi*F(t)) .* conj(F(t)) , 1) );
	evalP =@(y) 1/fc*real(F(t)'*y);
	
	p1 = evalP(F(x1)*a1(:));
	p2 = evalP(F(x2)*a2(:));
	p3 = evalP(F(x3)*a3(:));
	set(gca,'ylimspec','tight');
	plot(t,p1+p2+p3,'color',[0 0.4470 0.7410],'linewidth',2,'linestyle','-');
	xticks([]), yticks([]);
	box on

% 	m = floor(n/2);
%    O = ones([m+1,1]); O = padarray(O,m,'post'); % outer padding
%    D = fftshift(ifftn(abs(fftn(O)).^2));
% 	W = (n+1)*D(:).*D(:)';
% 	M1 = W.*(F(x1)*diag(a1)*F(x1)');
% 	M2 = W.*(F(x2)*diag(a2)*F(x2)');
% 	M3 = W.*(F(x3)*diag(a3)*F(x3)');
% 	p1 = evalP(M1);
% 	p2 = evalP(M2);
% 	p3 = evalP(M3);
% 	PP = (2*n+6)*(p1+p2+p3)/sum(p1+p2+p3);
% 	plot(t,PP,'r--','linewidth',2);


else
	scatter(x,y,ms,'.');
	xlim([0,1]), ylim([0 1])
end

lambda = [.3128 .6872];
lambda = lambda/sum(lambda);

%%


fc = 50;
n = [fc fc];
N = 2*n+1;
% frequency
[fY,fX] = meshgrid(0:fc,0:fc);
toepl = @(s) double(fY-fX==s);
 
% Fourier moments
F1 = @(x) exp(-2i*pi*(-fc:fc)'*x(:)');
F = @(x,y)exp( -2i*pi*(fX(:)*x(:)' + fY(:)*y(:)') );
% evaluate a polynomial associated to moment matrix on a grid
T = 256;
%t = linspace(0,1,m);
t = (0:T-1)/T;
[Y,X] = meshgrid(t,t); X = X(:); Y = Y(:);
evalP = @(Mi)real( sum( (Mi*F(X,Y)) .* conj(F(X,Y)) , 1) );
evalP = @(Mi)reshape(evalP(Mi),[T,T]);

% moment matrix
M1 = F1(x)*diag(a)*F1(x)';
M = F(x,y) * diag(a) * F(x,y)';
U = F(x,y) * sqrt(diag(a));

% % adding some noise
%sigma = 1e-2 * norm(M,'inf');
sigma = 0;
options.signed = 0;
[~,ids] = marginals(n,'colex',1);
mm = reshape(M(ids),N);
mm = ifftshift(mm);
%mm = mm + sigma*(randn(N)+1i*randn(N));
mm = mm + sigma*(randn(N));
mm(1) = real(mm(1));

mm = fftshift(mm);
momo =  genorder(n,'colex',1);
M_noisy = zeros(size(M));
for k1=-fc:fc
    for k2=-fc:fc
        Th = kron(toepl(k2), toepl(k1));
        M_noisy = M_noisy + mm(k1+fc+1,k2+fc+1) * Th;
    end
end
% M_noisy = conj(M_noisy+M_noisy')/2; % HACK?

% sigma = 1e-2;
% V  = rand(size(U))+1i*rand(size(U));
% V(1) = abs(V(1));
% U = U + sigma*V;

%

svals = svd(M);

% plot moment matrix svd
%clf, plot(log10(svals(1:150)),'.','markersize',15);
clf, plot(log10(svals),'.','markersize',15);
%xlabel('index')
ylabel('singular value (logarithm)')
%set(gca,'xtick',[])
box off
set(gca,'fontsize',21)

%mypath = '/home/pcatala/Workspace/Thesis/code/OT-SDP/results/moment_matrix_rank';
%saveas(gca,[mypath 'svd_sparse-measure-50-fc-25'],'epsc');

% Prony recovery

options.shift_mode  = 'kunis';
options.mode_debug = 0;
tol = 1e-1;
options.tol = tol;
options.jdiag = 'cardoso';
%options.lambda = lambda;
options.factorized = 0;

options.signed = 0;

%l = rand(2,1);
%options.lambda = 1;
%options.lambda = rand(2,1);

CA = (1-abs(a))*[.5 .5 1] + abs(a)*[1 .5 .5];

%profile on
%[supp,amp,modulus] = mvprony(M1,1,options);
[supp,amp,info] = mvprony(M_noisy,n,options);
amp = amp/sum(abs(amp));
%supp = prony2d(M,options);
%profile viewer

ell1x = supp(:,1); y1 = supp(:,2);
a1 = real( pinv(F(ell1x,y1)) * (F(x,y)*a) );

CA = (1-abs(a))*[.5 .5 1] + abs(a)*[1 .5 .5];
% m = rescale( abs(a1) ); m = m/max(m);
%m = (modulus(:,1)+modulus(:,2))/2;
m = info.modulus;
CM = (1-m)*[0 0 1] + m*[1 0 0];

%% display
scale=4e3;
clf, hold on;
%scatter(x1,y1,50,CM);
%plot(z,'r.','markersize',10);
%scatter(x,y,10*scale*a,CA,'filled','markerfacealpha',.3,'markeredgealpha',.3);
%scatter(ell1x,y1,scale*amp,CM,'filled');
scatter(x,y,1e5*a,CA,'o','linewidth',2.5);
scatter(1-ell1x,1-y1,1e3*amp,CM,'x','linewidth',3);
%scatter(ell1x,y1,100,CM,'filled');
%scatter(x1,y1,50,CM,'filled');
%plot([0 1],[0,-l(1)/l(2)],'k--')
%yy = mod(-l(1)/l(2) * t, 1);
%plot([0 1],[0 -l(1)/l(2)],'k--');
%plot(t,yy,'k--');
%plot([0 1],[2 0],'k--');
%plot(t,sqrt(t),'k--');
xlim([0,1]), ylim([0,1]);
xticks([]), yticks([]);
box on;
drawnow;

% for i=1:length(xalong)
%     ss = xalong{i};
%     x1 = ss(:,1); y1 = ss(:,2);
%     
%     clf, hold on;
%     scatter(x,y,'filled');
%     scatter(x1,y1,'filled');
%     xlim([0,1]),ylim([0 1]);
%     %pause(.05)
%     
%     drawnow;
% end


