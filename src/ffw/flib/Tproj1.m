function Tvals = Tproj(m,U)

r = size(U,2);
U = reshape(U,[m,r]);

PU = padarray(U,[m-1,0],'post');
Tvals = sum( ifft(abs(fft(PU)).^2), 2) ./ Dnumel1(m);
end
