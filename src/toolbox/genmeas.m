function measure = genmeas(d,measure)
% GENMEAS generate Radon measure

type = getoptions(measure,'type','discrete');

switch type
	case 'discrete'
		s = getoptions(measure,'sparsity',8);

		x = rand(s,d);
		a = rand(s,1);

	case 'gaussian'
		s   = getoptions(measure,'sparsity',500);
		mu  = getoptions(measure,'mean',0.5*ones(d,1));
		sig = getoptions(measure,'var',0.01);

		x = rand(s,d);
		a = rand(s,1);
end

a = a/norm(a,1);

measure.x = x;
measure.a = a;


end
