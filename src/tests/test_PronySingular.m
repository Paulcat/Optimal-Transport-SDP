% Test singular cases in random combination

curve = 'circle';
%curve = 'rand-curve';

n = 5000; % number of points
d = 2; % dimension

switch curve
	case 'circle'
		r = .2;
		z = .5+.5i + r*exp(2i*pi*linspace(0,1,n));
		
		x = real(z)';
		y = imag(z)';
		a = ones(n,1)/n;
		
	case 'rand-curve'
		% Randomize amplitude and phase.
		H = 20;
		rho = rand(1,H) .* logspace(-0.5,-2.5,H);
		phi = rand(1,H) .* 2*pi;
		
		% Accumulate r(t) over t=[0,2*pi]
		t = linspace(0,2*pi,n)';
		r = ones(size(t));
		for h=1:H
			r = r + rho(h)*sin(h*t+phi(h));
		end
		
		% Reconstruct z(t)
		x = .7 + .2*r .* cos(t);
		y = .3 + .2*r .* sin(t);
		a = ones(n,1)/n;
end

% random vector
la = rand(d,1);% + 1i*rand(d,1);
la = la/sum(la);
%la = [0.001;0.999];
%la = [.5;.5];

% "singular mapping"
arg = @(x) -angle(x)/2/pi;
f = @(x,y) mod(2*arg(sum(la.' .* exp(-2i*pi*[x(:),y(:)]),2)) - ...
	2*arg(la)' - [x(:),y(:)], 1);

Q = f(x,y);

%% prony recovery
fc = 50;
[fY,fX] = meshgrid(0:fc);
mm = exp(-2i*pi*(fX(:)*x(:)'+fY(:)*y(:)')) .* sqrt(a(:))';
options.factorized = 1;
options.jdiag = 'random';
options.lambda = la;
[u,s,info] = mvprony(mm,[fc,fc],options);

% barycenter
ba = la(1)*exp(-2i*pi*x) + la(2)*exp(-2i*pi*y);
id_small = find(abs(ba)<1e-3);
id_big = find(abs(ba)>.999);

%col1 = turbo(size(u,1));
col = turbo(size(x,1));

%% 
clf; hold on;
%subplot(1,2,1), hold on;
plot([0,1],[0,1],'k:','linewidth',2);
plot([0,.5],[.5,1],'k:','linewidth',2);
plot([0.5,1],[0,.5],'k:','linewidth',2);
xlim([0,1]),ylim([0,1]);
%
scatter(x,y,100,col,'.');
%plot(Q(:,1),Q(:,2),'--','color',[.7 .7 .7],'linewidth',2);
vv = [];
step = 40;
for ii = 0:2:n/step
	vv = [vv,(ii*step):((ii+1)*step)];
end
vv = vv(2:end);
%scatter(Q(:,1),Q(:,2),100,col,'.');
scatter(Q(vv,1),Q(vv,2),100,col(vv,:),'.');
%scatter(u(:,1),u(:,2),1e3*s,col1,'bx','linewidth',3);
%
%scatter(x(id_small),y(id_small),200,'r.');
xticks([]);
yticks([]);
box on;
axis tight
%scatter(x(id_big),y(id_big),200,'b.');

hack1 = plot([100 105],[100,105],'linewidth',3.5,'color',col(1000,:));
hack2 = plot([100 105],[100,105],'--','linewidth',3.5,'color',col(4500,:));
xlim([0,1]),ylim([0,1]);
%
legend([hack1,hack2],'Support $S$','Image $f_\lambda(S)$','location','southeast',...
	'interpreter','latex');
set(gca,'fontsize',30);

%%
% subplot(1,2,2), hold on;
clf, hold on;
eV = info.eVal;
col2 = parula(size(eV,1));
t = linspace(0,2*pi,100);
plot(cos(t),sin(t),'color',[.8 .8 .8],'linewidth',3);
plot([0,0],[-1,1],'--','color',[.8 .8 .8],'linewidth',3);
plot([-1,1],[0,0],'--','color',[.8 .8 .8],'linewidth',3);
%scatter(real(ba),imag(ba),100,col,'.');
scatter(real(ba),imag(ba),50,[.6 .6 .6],'.');
sc = scatter(real(eV),imag(eV),600,col2,'.');
xticks([]);
yticks([]);
axis tight;
box on;

%hack1 = plot([100 105],[100,105],'linewidth',3.5,'color',col(1000,:));
hack1 = scatter(100,100,100,col(1000,:),'filled');
xlim([-1,1]),ylim([-1,1]);
%legend(hack1,'Image $\psi_\lambda(S)$','location','southeast',...
%	'interpreter','latex');
legend(hack1,'Eigenvalues of $X_\lambda$','location','southeast',...
	'interpreter','latex');
set(gca,'fontsize',30);

%% reconstruction

clf, hold on;
plot([0,1],[0,1],'k:','linewidth',2);
plot([0,.5],[.5,1],'k:','linewidth',2);
plot([0.5,1],[0,.5],'k:','linewidth',2);
scatter(x,y,100,[.9 .9 .9],'.');
scatter(Q(vv,1),Q(vv,2),100,[.9 .9 .9],'.');
scatter(u(:,1),u(:,2),1e4*s,col2,'.','linewidth',3);
xlim([0,1]);
ylim([0,1]);
xticks([]);
yticks([]);
box on;

hack1 = plot([100,105],[100 105],'linewidth',3.5,'color',[.9 .9 .9]);
hack2 = scatter(100,100,100,col(1000,:),'filled');
legend([hack1,hack2],'True support $S$','Reconstructed support','location','southeast',...
	'interpreter','latex');
set(gca,'fontsize',30);

%% dynamic drawing
% figure(1), clf, hold on;
% plot([0,1],[0,1],'k','linewidth',2);
% xlim([0,1]),ylim([0,1]);
% 
% figure(2), clf;
% ax2H = axes('nextplot','add'); % = hold on
% plot(cos(t),sin(t),'k','linewidth',2);
% plot([0,0],[-1,1],'k--','linewidth',2);
% plot([-1,1],[0,0],'k--','linewidth',2);
% xlim([-1.5,1.5]),ylim([-1.5,1.5]);
% 
% col = turbo(size(x,1));
% for i=1:size(x,1)
% 	figure(1);
% 	scatter(x(i),y(i),100,col(i,:),'.');
% 	scatter(Q(i,1),Q(i,2),100,col(i,:),'.');
% 	drawnow;
% 	
% 	figure(2);
% 	delete([sc1,sc2]);
% 	sc1 = scatter(cos(2*pi*x(i)),-sin(2*pi*x(i)),100,'b.');
% 	sc2 = scatter(cos(2*pi*y(i)),-sin(2*pi*y(i)),100,'r.');
% 	z0 = la(1)*exp(-2i*pi*x(i)) + la(2)*exp(-2i*pi*y(i));
% 	sc0 = scatter(real(z0),imag(z0),100,'g.');
% 	%
% end