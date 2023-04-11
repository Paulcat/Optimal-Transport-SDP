% example OT between 2D discrete measures

clear all

% *** set up problem ***
% **********************

% create marginals
[sA,sB] = deal(2,2);
% 1st
xA = rand(sA,2);
a  = rand(sA,1);
a  = a/norm(a,1);
% 2nd
xB = rand(sB,2);
b  = rand(sB,1);
b  = b/norm(b,1);


%% moments of marginals
[nA1,nA2,nB1,nB2] = deal(4,4,4,4);
nn  = [nA1,nA2,nB1,nB2];
KA  = genorder([nA1,nA2],'colex',0); % same as [fY1,fX2] = meshgrid(-n1:n1,-n2:n2)
KB  = genorder([nB1,nB2],'colex',0);
KAp = genorder([nA1,nA2],'colex',1);
KBp = genorder([nB1,nB2],'colex',1);
% 1st marginal (A)
cA = exp(-2i*pi * (KA (:,1)*xA(:,1)' + KA (:,2)*xA(:,2)')) * a(:);
FA = exp(-2i*pi * (KAp(:,1)*xA(:,1)' + KAp(:,2)*xA(:,2)'));
UA = FA .* sqrt(a');
% 2nd marginal (B)
cB = exp(-2i*pi * (KB (:,1)*xB(:,1)' + KB (:,2)*xB(:,2)')) * b(:);
FB = exp(-2i*pi * (KBp(:,1)*xB(:,1)' + KBp(:,2)*xB(:,2)'));
UB = FB .* sqrt(b');
% shape and shift (order 0 comes first)
cA = reshape(cA,[2*nA1+1, 2*nA2+1]);
cA = ifftshift(cA);
cB = reshape(cB,[2*nB1+1, 2*nB2+1]);
cB = ifftshift(cB);


% polynomial display
[NA1,NA2,NB1,NB2] = deal(512,512,512,512);
[tA1,tA2] = ndgrid((0:NA1-1)'/NA1, (0:NA2-1)'/NA2);
[tB1,tB2] = ndgrid((0:NB1-1)'/NB1, (0:NB2-1)'/NB2);
%
evalA = @(p) reshape(real( exp(-2i*pi * (KA(:,1)*tA1(:)' + KA(:,2)*tA2(:)'))'*p(:) ), [NA1,NA2]);
evalB = @(p) reshape(real( exp(-2i*pi * (KB(:,1)*tB1(:)' + KB(:,2)*tB2(:)'))'*p(:) ), [NB1,NB2]);

% hyper parameters for ffw
Cl = (norm(evalA(cA),'inf') + norm(evalB(cB),'inf'))/2;
Cr = prod(nn) * (norm(cA,'inf') + norm(cB,'inf'))/2; %??
%
la  = 1e-5*Cl;
rho = 1e-5*Cr;



% transport cost (in Fourier domain)
cost = zeros(2*nn+1);
cost(nA1+1, nA2+1, nB1+1, nB2+1) = 1;
cost(nA1+2, nA2+1, nB1  , nB2+1) = -1/4;
cost(nA1+1, nA2+2, nB1+1, nB2  ) = -1/4;
cost(nA1  , nA2+1, nB1+2, nB2+1) = -1/4;
cost(nA1+1, nA2  , nB1+1, nB2+2) = -1/4;
% shift
cost = ifftshift(cost);


% variable size and functionals
mm 	 = nn+1;
[f,f0] = ot_fobj (mm,cost,cA,cB,la);
g 		 = ot_fgrad(mm,cost,cA,cB,f0,la);
Tpen = ffw_Tpen(mm);


% **** load problem ****
% **********************
problem = struct;
%
problem.name 	= 'optimal-transport';
problem.vardim = mm;
problem.fobj 	= f;
problem.f0 		= f0;
problem.grad 	= g;
problem.cflag  = 'none';
problem.hyper 	= la;
problem.ls 		= ot_lscoeffs(mm,cost,cA,cB,f0,la);
% **********************
% **********************



% **** load solver options ****
% *****************************
options = struct;
%
options.tol 			= 1e-5;
options.maxiter 		= 20;
options.bfgsToolbox  = 'minfunc';
options.bfgsProgTol 	= 1e-16;
options.bfgsMaxIter 	= 500;
options.lmoTol 		= 1e-10;
options.lmoMaxIter 	= 1e3;
options.rho 			= rho;
options.display 		= 'on';
% ******************************
% ******************************



% *** cvx ground truth ***
% ************************
% periodic cost
XA = reshape(xA,[sA,1,2]);
XB = reshape(xB,[1,sB,2]);
C 	= sum( sin(pi*(XA-XB)).^2, 3);

% cvx solver
cvx_begin
variable P0(sA,sB) nonnegative
%
sum(P0,1) == b.';
sum(P0,2) == a;
%
minimize( C(:)'*P0(:) );
cvx_end
val_cvx = C(:)'*P0(:);


% recovery
[I,J] = find(P0 > 1e-4 * max(abs(P0(:))));
x0 = [xA(I,:),xB(J,:)];
a0 = P0(sub2ind([sA,sB],I,J));
%
Freq = genorder(nn,'colex',1);
Freq = reshape(Freq,[prod(nn+1),1,4]);
Pos  = reshape(x0,[1,size(x0,1),4]);
F0 = exp(-2i*pi * sum( Freq .* Pos, 3));
U0 = F0 .* sqrt(a0');
T0 = Tprojn(mm,U0);
%
c0A = exp(-2i*pi * (KA(:,1)*x0(:,1)' + KA(:,2)*x0(:,2)')) * a0(:);
c0A = ifftshift(reshape(c0A,[2*nA1+1,2*nA2+1]));
c0B = exp(-2i*pi * (KB(:,1)*x0(:,1)' + KB(:,2)*x0(:,2)')) * a0(:);
c0B = ifftshift(reshape(c0B,[2*nB1+1,2*nB2+1]));

% ***********************
% ***********************



% *** OUR SOLVER ***
% ******************
[U,info] = FFW(problem,options);

% prony extraction
options_prony.factorized = 1;
options_prony.jdiag = 'cardoso';
[x,amp] = mvprony(U,nn,options_prony);


% objective values
val_ffw = info.E(end);
val_gt_in_penalized   = f(T0) + 1/rho*Tpen(U0,T0);
%
CC = sum( sin(pi*(x(:,1:2)-x(:,3:4))).^2, 2);
val_ffw_in_constrained = sum(CC(:) .* amp(:));

%% marginals
mA = exp(-2i*pi * (KA(:,1)*x(:,1)' + KA(:,2)*x(:,2)')) * amp(:);
mA = ifftshift(reshape(mA,[2*nA1+1,2*nA2+1]));
mB = exp(-2i*pi * (KB(:,1)*x(:,3)' + KB(:,2)*x(:,4)')) * amp(:);
mB = ifftshift(reshape(mB,[2*nB1+1,2*nB2+1]));
err_mA = norm(mA-cA,'fro')/norm(cA,'fro');
err_mB = norm(mB-cB,'fro')/norm(cB,'fro');

% toeplitz
err_T = sqrt(2*Tpen(U,Tprojn(mm,U))) / norm(U'*U,'fro');


fprintf('INFOS\n');
fprintf('%35s : %-10.7f\n', 'FFW final value', val_ffw);
fprintf('%35s : %-10.7f\n', 'GT (CVX) solution in penalized obj', val_gt_in_penalized);
fprintf('%35s : %-10.7f\n', 'FFW solution in constrained obj', val_ffw_in_constrained);
fprintf('%35s : %-10.7f\n', 'CVX final value', val_cvx);
fprintf('%35s : %-10.7f\n', '1st marginal constraint violation', err_mA);
fprintf('%35s : %-10.7f\n', '2nd marginal constraint violation', err_mB);
fprintf('%35s : %-10.7f\n', 'Toeplitz constraint violation', err_T);


% display
clf, hold on;
scatter(x0(:,1),x0(:,2),50,'r','filled');
scatter(x0(:,3),x0(:,4),50,'b','filled');
for i=1:size(x,1)
	plot([x0(i,1),x0(i,3)], [x0(i,2),x0(i,4)], 'k', 'linewidth', 5*a0(i));
end

scatter(x(:,1),x(:,2),50,'g','filled');
scatter(x(:,3),x(:,4),50,'m','filled');
for i=1:size(x,1)
	plot([x(i,1),x(i,3)], [x(i,2),x(i,4)], 'c', 'linewidth', 5*a0(i));
end
