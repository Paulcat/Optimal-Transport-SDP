function Tx = Tprodn(m,Tvals,x)
% TPRODN n-level Toeplitz product


d = length(m);
r = size(x,2);
x = reshape(x,[m,r]);

Px = padarray(x,[m-1,0],'post');
Cx = ifftd(d, fftd(d,Tvals) .* fftd(d,Px) ); % broadcasting

% cropping in arbitrary dimension requires cells...
crop = cellfun(@(j) 1:j, num2cell(m), 'un', 0);
Tx = Cx(crop{:},:);
Tx = reshape(Tx,prod(m),r);
