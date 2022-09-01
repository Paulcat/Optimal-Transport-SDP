function [Tpen,Tpen_g,Tproj,Tprod,Dnumel] = ffw_Tpen(m)
% FFW_TPEN Toeplitz penalization term

d = numel(m);

switch d
	case 1
		Tproj  = @(U) 	 Tproj1(m,U);
		Tprod  = @(T,v) Tprod1(m,T,v);
		Dnumel = Dnumel1(m);

	case 2
		Tproj  = @(U)   Tproj2(m,U);
		Tprod  = @(T,v) Tprod2(m,T,v);
		Dnumel = Dnumel2(m);

	case 4
		Tproj  = @(U)   Tproj4(m,U);
		Tprod  = @(T,v) Tprod4(m,T,v);
		Dnumel = Dnumel4(m);
end

% Toeplitz penalized term
normT_sq = @(T) 		sum( Dnumel .* abs(T).^2, 1:d); % Toeplitz norm
Tpen 		= @(U,TU) 	1/2 * (norm(U'*U,'fro')^2 - normT_sq(TU)); % penalization
Tpen_g 	= @(U,TU,h) U*(U'*h) - Tprod(TU,h); % gradient

end
