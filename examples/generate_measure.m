function [X,A] = generate_measure(type)

switch type
	case 'sparse'
		s = 10;

		x = rand(s,1);
		y = rand(s,1);
		X = [x(:), y(:)];


	case 'curve-circle'
		s = 1000;
		t = linspace(0,1,s);
		%t = rand(1,s);

		r0 = .3;
		[c0x,c0y] = deal(.5,.5);

		x = c0x + r0 .* cos(2*pi*t);
		y = c0y + r0 .* sin(2*pi*t);
		X = [x(:), y(:)];


	case 'curve-treffle'
		s = 1000;
		t = linspace(0,1,s);

		% center and radius
		[c0x,c0y] = deal(.5,.5);
		r0 = .35;

		be = 12;
		z = .15 * (1.3 + cos(be*pi*t)) .* exp(2i*pi*t);
		z = r0 * z/max(abs(z));
		z = z + c0x + c0y*1i;

		x = real(z);
		y = imag(z);
		X = [x(:), y(:)];


	case 'curve-lissajou1'
		s = 1000;
		t = linspace(0,1,s);

		% center and ~radius
		[c0x,c0y] = deal(.5,.1);
		r0 = 3/4;

		x = c0x + r0 * cos(pi*t).^3 .* sin(pi*t);
		y = c0y + r0 * cos(pi*t).^2;
		X = [x(:), y(:)];


	case 'curve-rand'
		s = 1000;
		t = linspace(0,2*pi,s);

		% center and ~radius
		[c0x,c0y] = deal(.5,.5);
		r0 = .2;

		% random amplitude and phase
		H = 20;
		rho = rand(1,H) .* logspace(-.5,-2.5,H);
		phi = rand(1,H) .* 2*pi;

		% accumulate r(t) over t=[0,2*pi]
		r = ones(size(t));
		for h=1:H
			r = r + rho(h) * sin(h*t + phi(h));
		end
		r = r0 * r;

		x = c0x + r .* cos(t);
		y = c0y + r .* sin(t);
		X = [x(:), y(:)];


	case 'curve-rand-3d'


	case 'area-disc'
		s = 10000;
		t = linspace(0,1,s);

		% center and radius
		[c0x,c0y] = deal(.5,.5);
		r0 = .2;

		% uniform sampling?
		th = 2*pi*rand(s,1)
		r  = rand(s,1);

		% poisson sampling?

		% regular sampling?

		x = c0x + r0 * sqrt(r) .* cos(th);
		y = c0y + r0 * sqrt(r) .* sin(th);
		X = [x(:), y(:)];


	case 'area-rand'
		s = 10000; % lower bound on the number of points

		% center and ~radius
		[c0x,c0y] = deal(.5,.5);
		r0 = .2;

		% randomize amplitude and phase
		H = 20;
		rho = rand(1,H) .* logspace(-.5,-2.5,H);
		phi = rand(1,H) .* 2*pi;

		% fill area radially
		nrad = 100;
		al	  = ceil(2*s/nrad/(nrad+1));
		s    = al * nrad * (nrad+1) / 2; % actual number of points
		%
		x = [];
		y = [];
		for i=1:nrad
			t = linspace(0,2*pi,al*i);

			% accumulate r(t) over t=[0,2*pi]
			r = ones(size(t));
			for h=1:H
				r = r + rho(h) * sin(h*t + phi(h));
			end

			% fill up area
			x = [x, c0x + i/nrad * r0 * (r + .15*rand(size(r))) .* cos(t + .15*rand(size(t)))];
			y = [y, c0y + i/nrad * r0 * (r + .15*rand(size(r))) .* sin(t + .15*rand(size(t)))];
		end
		X = [x(:), y(:)];


	case 'sphere'
		s = 10000;

		% center and radius
		[c0x,c0y,c0z] = deal(.5,.5,.5);
		r0 = .3;

		th  = acos(1-2*rand(s,1));
		phi = 2*pi * rand(s,1);

		x = c0x + r0 * sin(th).*cos(phi);
		y = c0y + r0 * sin(th).*sin(phi);
		z = c0z + r0 * cos(th);
		X = [x(:), y(:), z(:)];

end

A = ones(s,1)/s; % uniform weights...

end
