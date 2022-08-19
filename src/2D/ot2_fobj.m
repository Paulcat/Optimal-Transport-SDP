function F = ot2_fobj(m,cost,u1,u2,la,rho)
%FOBJ OT objective for FFW (2D)

debug = 0;

% scaling
%f0 = 4*la  / (norm(u1(:),'fro')^2 + norm(u2(:),'fro')^2);
f0 = 1;

% helper
d = length(m);
flat = @(x)x(:);
normT2 = @(T) sum( Dnumel4(m) .* abs(T).^2, 1:d);


% objective part: terms depending on TVALS
FT = @(T) real(cost(:)'*T(:)) + ...
    1/4/la * ( norm(flat(T(:,:,1,1))-u1(:),'fro')^2 + ...
               norm(flat(T(1,1,:,:))-u2(:),'fro')^2 ) + ...
    -1/2/rho * normT2(T);

% objective
F = @(U) f0 * ( FT(Tproj4(m,U)) + 1/2/rho * norm(U'*U,'fro')^2);

if debug
    disp('debug mode!!');
    FT = @(T) 1/2*normT2(T);
    F  = @(U) FT(Tproj4(m,U));
end

end

