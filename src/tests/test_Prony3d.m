% Test Prony 3D

% generate a curve
name = 'rand-curve';
%name = 'lissajou';
%name = 'sphere';

switch name
    case 'rand-curve'
        % TODO
        n = 1000;
        
        % Randomize amplitude and phase
        H = 20;
        rho = rand(1,H) .* logspace(-0.5,-2.5,H);
        phi = rand(1,H) .* 2*pi;
        
        % Accumulate r(t) over t=[0,2*pi]
        t = linspace(0,2*pi,n);
        p = linspace(0,pi,n);
        r = ones(size(t));
        for h=1:H
            r = r + rho(h)*sin(h*t+phi(h));
        end
        
        % Reconstruct z(t)
        x = .5 + .3*r .* sin(t) .* cos(p);
        y = .5 + .3*r .* sin(t) .* sin(p);
        z = .5 + .3*r .* cos(t);
        
        a = ones(n,1)/n;
        
    case 'lissajou'
        n = 1000;
        t = linspace(0,2*pi,n);
        
        x = .5 + .1*4*sin(t+pi/2);
        y = .5 + .1*2*sin(3*t);
        z = .5 + .1*5*sin(2*t);
        
        %a = rand(n,1)/n;
		  a = ones(n,1)/n;
        
    case 'sphere'
        n = 4000;
        
        r = .3;
        theta = acos(1-2*rand(n,1));
        phi = 2*pi * rand(n,1);
        
        x = .5 + r*sin(theta).*cos(phi);
        y = .5 + r*sin(theta).*sin(phi);
        z = .5 + r*cos(theta);
        
        a = ones(n,1)/n;
end

ms = 50;

clf;
scatter3(x,y,z,ms,'.');
xlim([0,1]),ylim([0,1]),zlim([0,1]);

%%

fc = 15;
n = [fc fc fc];
[fX,fY,fZ] = ndgrid(0:fc);

% Fourier moments
F = @(x,y,z) exp(-2i*pi*(fX(:)*x(:)'+fY(:)*y(:)'+fZ(:)*z(:)'));

% moment matrix
M = F(x,y,z)*diag(a)*F(x,y,z)';
U = F(x,y,z) * sqrt(diag(a));

%% Prony recovery

options.shift_mode  = 'kunis';
options.mode_debug = 0;
tol = 1e-3;
options.tol = tol;
options.jdiag = 'cardoso';
options.lambda = rand(1,3);
options.factorized = 1;
options.signed = 0;

[supp,amp,info] = mvprony(U,n,options);
x1 = supp(:,1); y1 = supp(:,2); z1 = supp(:,3);
a1 = real( pinv(F(x1,y1,z1)) * (F(x,y,z)*a) );
amp = amp/sum(abs(amp));


%% display
CA = (1-abs(a))*[.5 .5 1] + abs(a)*[1 .5 .5];
m = info.modulus;
CM = (1-m)*[0 0 1] + m*[1 0 0];

[X,Y,Z] = sphere(100);
X = .5 + .29*X;
Y = .5 + .29*Y;
Z = .5 + .29*Z;

scale = 80e2;
clf;
%scatter3(x,y,z,10,CA,'filled','markerfacealpha',.3,'markeredgealpha',.3);
surf(X,Y,Z,'linestyle','none','facecolor',[.7 .7 1],'facealpha',1);
hold on;
scatter3(x1,y1,z1,scale*amp,CM,'filled');
%scatter3(x1,y1,z1,50,CM,'filled');
xlim([0,1]), ylim([0,1]), zlim([0,1]);
xticks([]), yticks([]), zticks([]);
box on;
drawnow;