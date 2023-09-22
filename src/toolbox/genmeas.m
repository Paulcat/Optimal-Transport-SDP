function [x,a] = genmeas(type,d,s,varargin)
% GENMEAS generate Radon measure

switch type
   case 'discrete'
      %s = getoptions(measure,'sparsity',8);
      
      x = rand(s,d);
      a = rand(s,1); % TODO: complex amplitudes?
      
   case 'gaussian'
      %s   = getoptions(measure,'sparsity',500);
      mu  = getoptions(measure,'mean',0.5*ones(d,1));
      sig = getoptions(measure,'var',0.01);
      
      x = rand(s,d);
      a = rand(s,1);
end

a = a/norm(a,1);

%measure.x = x;
%measure.a = a;


end
