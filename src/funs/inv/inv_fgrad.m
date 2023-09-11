function G = inv_fgrad(m,L,f0,la)
%INV_FGRAD Gradient for invariant measure problem
%
%   See also INV_FOBJ.

% A = ifftshift(L'*L);
d = numel(m);
A = L'*L;
A = reshape(A,[2*m-1, 2*m-1]);
A = permute(A,[d+1:2*d,1:d]);
prodA = @(T) shiftdim(sum(T.*A,1:d)); % remove all leading singleton dimensions
warning('TODO: check everything');

N = 1./Dnumeln(m);

G = @(T,h) 1/f0 * Tprodn(m,N.*( 1/la*prodA(T) + T), h);

end
