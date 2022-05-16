% Test multivariate Prony on a few 2D examples

type = 'sparse';
%type = 'curve-circle';
%type = 'curve-treffle';
%type = 'curve-lissajou1';
%type = 'curve-rand';
%type = 'area-disc';
%type = 'area-rand';



% helpers
d  = 2; % dimension
n 	= 10;
nn = [n n];
fX = compute_grid(nvec,'spectral');
fX	= reshape(fX,prod(nn+1),d);
Fn = @(x) exp(-2i*pi*(fX(:,1)*x(:,1)' + fX(:,2)*x(:,2)') );


% compute moment matrix
[x,a] = generate_measure(type);
U = Fn(x) .* sqrt(a');
M = U*U';


% recovery
options.factorized = 1;
[xrec,arec] = mvprony(U,nn,d,options);


clf, hold on;
scatter(x(:,1),x(:,2),50,'.');
scatter(xrec(:,1),xrec(:,2),50,'x','linewidth',2);
xlim([0,1]), ylim([0,1]);
