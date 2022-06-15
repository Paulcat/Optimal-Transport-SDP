function checkgradient(f,g,x)
%CHECKGRAD Numerical validation for gradient
%   Detailed explanation goes here

ndir = 100; % nb of directional derivatives
e    = 1e-4;

df      = zeros(1,ndir);
df_u = zeros(1,ndir); % user gradient
for i=1:ndir
    dir = rand(size(x)) ;%+ 1i*rand(size(x));
    dir = dir / norm(dir,'fro');
    
    df(i)   = f(x + e*dir) - f(x);
    df_u(i) = trace( g(x)' * dir );
end
df = df / e;


clf;
subplot(1,2,1), hold on;
plot(1:ndir, real(df),   '.');
plot(1:ndir, real(df_u), 'o');
%
subplot(1,2,2), hold on;
plot(1:ndir, imag(df),   '.');
plot(1:ndir, imag(df_u), 'o');

%df_u
%df ./ df_u


end
