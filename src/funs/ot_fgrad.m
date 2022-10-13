function G = ot_fgrad(m,cost,u1,u2,f0,la)
%OT_FGRAD Gradient of mmd-penalized ot (in arbitrary dimension)
%
%   See also OT_FOBJ

debug = 0;


% dimensions
d  = length(m); % dimension of transport plan
dm = d/2; % dimension of marginals

% safety net(s)
if ~prod(size(u1,1:dm)==2*m(1:dm)-1) || ~prod(size(u2,1:dm)==2*m(dm+1:d)-1)
	error('moment vector incorrectly shaped');
end



% gradient 'subroutines'
% I1(:,:,...,1,1,...) = 1 for arbitrary dimension
k1 = cell(1,d);
k1(1:dm) = {':'}; k1(dm+1:d) = {1};
I1 = zeros(2*m-1); I1(k1{:}) = 1;

% I2(1,1,...,:,:,...) = 1 for arbitrary dimension
k2 = cell(1,d);
k2(1:dm) = {1}; k2(dm+1:d) = {':'};
I2 = zeros(2*m-1); I2(k2{:}) = 1;

I 		= I1+I2;
u2t 	= reshape(u2,[ones(1,dm),size(u2)]);
Y 		= u1.*I1 + u2t.*I2;




% wrapper for gradient
N = 1./Dnumeln(m);
G = @(T,h) 1/f0 * ( Tprodn(m,N.*cost,h) + 1/2/la * Tprodn(m,N.*(T.*I-Y),h) );



if debug
    disp('debug mode!!');
    GT = @(T,h) 2*Tprodn(m,N.*(T.*I-Y),h);
    G  = @(U,h) Tprod2(m,Tproj2(m,U),h);
end

end

