function Tx = Tprod(m,Tvals,x)

r = size(x,2);
x = reshape(x,[m,r]);

Px = padarray(x,[m-1,0],'post');
Cx = ifft( fft(Tvals) .* fft(Px) );

Tx = Cx(1:m,:);
end