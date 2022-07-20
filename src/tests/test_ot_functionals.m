% Test for OT ffw functionals...
clear all

s = 500;
t = linspace(0,1,s);

[n1,n2] = deal(6,6);
nn 	  = [n1,n2];


% first measure
mu1 = .001 + .3 * exp(-(t-.3).^2 ./ (.02)^2 ) + .5 * exp(-(t-.6).^2 ./ (.08)^2 );
mu1 = mu1/norm(mu1,1);
%
%M1 = exp(-2i*pi*(0:n)'*t(:)') .* sqrt(mu1);
%M1 = M1*M1'; % moment matrix
c1 = exp(-2i*pi*(-n1:n1)'*t(:)') * mu1(:); % moment vector
c1 = ifftshift(c1);


% second measure
mu2 = .001 + .4 * exp(-(t-.5).^2 ./ (.04)^2 );
mu2 = mu2/norm(mu2,1);
%
%M2 = exp(-2i*pi*(0:n)'*t(:)') .* sqrt(mu2);
%M2 = M2*M2'; % moment matrix
c2 = exp(-2i*pi*(-n2:n2)'*t(:)') * mu2(:); % moment vector
c2 = ifftshift(c2);

% polynomials
[N1,N2] = deal(512,512);
[t1,t2] = deal( (0:N1-1)'/N1, (0:N2-1)'/N2 );
NN 	  = [N1,N2];
%
eval1 = @(p) real( exp(-2i*pi*(-n1:n1)'*t1')' * p(:) );
eval2 = @(p) real( exp(-2i*pi*(-n2:n2)'*t2')' * p(:) );


% cost
h = zeros(2*nn+1);
h(n1:n1+2, n2:n2+2) = [0 0 -1/4; 0 1/2 0; -1/4 0 0];
h = ifftshift(h);

%% test gradient

% parameters
[la,rho] = deal(rand,rand);


[f,f0] = ot1_fobj (nn+1,h,c1,c2,la,rho); % TODO: why do I force n+1 here?
g		 = ot1_fgrad(nn+1,h,c1,c2,f0,la,rho);
g		 = @(u) 2*g(u,u);
U = rand(prod(nn+1),6) + 1i*rand(prod(nn+1),6);
checkgradient(f,g,U);

%clf, hold on;
%area(t,mu1);
%area(t,mu2);


%% test scalings

if ~exist('v','var')
	v = rand(prod(nn+1),1) + 1i*rand(prod(nn+1),1);
end

Cl = (norm(eval1(c1),'inf') + norm(eval2(c2),'inf'))/2;
Cr = prod(nn) * (norm(c1,'inf') + norm(c2,'inf'))/2;
%
la_list  = Cl * [1e-10,1e-6,1e-3,1e-1,1e1,1e3,1e6,1e10];
rho_list = Cr * [1e-10,1e-6,1e-3,1e-1,1e1,1e3,1e6,1e10];

testf = zeros(length(la_list),length(rho_list),2);
testg = zeros(length(la_list),length(rho_list),2);
%
testf_r = zeros(length(la_list),length(rho_list),2);
testg_r = zeros(length(la_list),length(rho_list),2);
%
testf1 = zeros(length(la_list),length(rho_list),2);
testg1 = zeros(length(la_list),length(rho_list),2);
%
testf2 = zeros(length(la_list),length(rho_list),2);
testg2 = zeros(length(la_list),length(rho_list),2);
%
testf3 = zeros(length(la_list),length(rho_list),2);
testg3 = zeros(length(la_list),length(rho_list),2);
%
f0s	= zeros(length(la_list),length(rho_list));
%
U = rand(prod(nn+1),6) + 1i*rand(prod(nn+1),6);

%la = Cl*1e-3;
%rho = Cr*1e1;
for i = 1:length(la_list)
	fprintf('%d,%',i); pause(1e-10); %hack?
	la = la_list(i);
	
	for j = 1:length(rho_list)
		rho = rho_list(j);
		
		[f,f0] = ot1_fobj (nn+1,h,c1,c2,la,rho);
		[g,gt] = ot1_fgrad(nn+1,h,c1,c2,f0,la,rho);
		
		% functions parts
		D  = Dnumel2(nn+1);
		f1 = @(T) 1/4/f0/la*(norm(T(:,1)-c1,'fro')^2 + norm(T(1,:).'-c2,'fro')^2);
		f2 = @(U,T) 1/2/f0/rho*(norm(U'*U,'fro')^2 - sum(D.*abs(T).^2, 1:2));
		f3 = @(T) 1/f0*real(trace(h'*T));
		%
		I1 = zeros(2*nn+1); I1(:,1) = 1;
		I2 = zeros(2*nn+1); I2(1,:) = 1;
		I  = I1+I2;
		g1 = @(T,w) 1/2/f0/la * Tprod2(nn+1,1./D.*(T.*I-c1.*I1-c2.'.*I2),w);
		g2 = @(U,T,w) 1/f0/rho * ( U*(U'*w) - Tprod2(nn+1,T,w) );
		g3 = @(w) 1/f0 * Tprod2(nn+1,1./D.*h,w);
		
		% at solution
		options.maxiter = 500;
		options.display = 'off';
		U0 = ffw_bfgs(U,f,@(U)2*g(U,U),options);
		%
		T  = Tproj2(nn+1,U);
		T0 = Tproj2(nn+1,U0);
		
		% sanity check
		%s1 = f(U) - f1(T)-f2(U,T)-f3(T)
		%s2 = norm(g(U,v) - g1(T,v)-g2(U,T,v)-g3(v))
		
		testf(i,j,:) = [f(U),f(U0)];
		testg(i,j,:) = [norm(g(U,v)), norm(g(U0,v))];
		%
		testf_r(i,j,:) = rho/max(la,rho) * [f(U),f(U0)];
		testg_r(i,j,:) = rho/max(la,rho) * [norm(g(U,v)), norm(g(U0,v))];
		%
		testf1(i,j,:) = [f1(T), f1(T0)];
		testf2(i,j,:) = [f2(U,T), f2(U0,T0)];
		testf3(i,j,:) = [f3(T), f3(T0)];
		%
		testg1(i,j,:) = [norm(g1(T,v)),norm(g1(T0,v))];
		testg2(i,j,:) = [norm(g2(U,T,v)),norm(g2(U,T0,v))];
		testg3(i,j,:) = [norm(g3(v)),norm(g3(v))];
		%
		f0s(i,j)	= f0;
	end
end
fprintf('\n');


% display 1d 
clf;
subplot(2,2,1);
loglog(la_list,testf1(:,3,1),'c:','linewidth',2);
hold on;
loglog(la_list,testf2(:,3,1),'c-.','linewidth',2);
loglog(la_list,testf3(:,3,1),'c--','linewidth',2);
%
loglog(la_list,testf(:,3,1),'b','linewidth',2);
loglog(la_list,testf_r(:,3,1),'b-o','linewidth',2);
%
%loglog(la_list,f0s,'color',[.5 .5 .5],'linewidth',2);
%loglog(la_list,1./la_list,'--','color',[.5 .5 .5],'linewidth',2);
legend('terme en \lambda','terme en \rho','terme coût','total','rescaled total','location','southeast');
title('FOBJ: at random point');

subplot(2,2,2);
loglog(la_list,testg1(:,3,1),'m:','linewidth',2);
hold on;
loglog(la_list,testg2(:,3,1),'m-.','linewidth',2);
loglog(la_list,testg3(:,3,1),'m--','linewidth',2);
%
loglog(la_list,testg(:,3,1),'r','linewidth',2);
loglog(la_list,testg_r(:,3,1),'r-o','linewidth',2);

legend('terme en \lambda','terme en \rho','terme coût','total','rescaled total','location','southeast')
title('Norme du GRAD: at random point');

subplot(2,2,3);
loglog(la_list,testf1(:,3,2),'c:','linewidth',2);
hold on;
loglog(la_list,testf2(:,3,2),'c-.','linewidth',2);
loglog(la_list,testf3(:,3,2),'c--','linewidth',2);
%
loglog(la_list,testf(:,3,2),'b','linewidth',2);
loglog(la_list,testf_r(:,3,2),'b-o','linewidth',2);
%
%loglog(la_list,f0s,'color',[.5 .5 .5],'linewidth',2);
%loglog(la_list,1./la_list,'--','color',[.5 .5 .5],'linewidth',2);

legend('terme en \lambda','terme en \rho','terme coût','total','rescaled total','location','southeast');
title('FOBJ: after BFGS minimization');

subplot(2,2,4);
loglog(la_list,testg1(:,3,2),'m:','linewidth',2);
hold on;
loglog(la_list,testg2(:,3,2),'m-.','linewidth',2);
loglog(la_list,testg3(:,3,2),'m--','linewidth',2);
%
loglog(la_list,testg(:,3,2),'r','linewidth',2);
loglog(la_list,testg_r(:,3,2),'r-o','linewidth',2);
legend('terme en \lambda','terme en \rho','terme coût','total',' rescaled total','location','southeast')
title('Norme du GRAD: after BFGS minimization')

% display 2d
clf;
[R,L] = ndgrid(la_list,rho_list);

fs = 6;

subplot(2,2,1);
surf(L,R,log(testg(:,:,1)),'linestyle','none');
set(gca,'xscale','log');
set(gca,'yscale','log');
colorbar
%
xlabel('\lambda')
ylabel('\rho');
set(gca,'fontsize',fs);

subplot(2,2,2);
surf(L,R,log(testg(:,:,2)),'linestyle','none');
set(gca,'xscale','log');
set(gca,'yscale','log');
colorbar
%
xlabel('\lambda')
ylabel('\rho');
set(gca,'fontsize',fs);

subplot(2,2,3);
surf(L,R,log(testg_r(:,:,1)),'linestyle','none');
set(gca,'xscale','log');
set(gca,'yscale','log');
colorbar
%
xlabel('\lambda')
ylabel('\rho');
set(gca,'fontsize',fs);

subplot(2,2,4);
surf(L,R,log(testg_r(:,:,2)),'linestyle','none');
set(gca,'xscale','log');
set(gca,'yscale','log');
colorbar
%
xlabel('\lambda')
ylabel('\rho');
set(gca,'fontsize',fs);