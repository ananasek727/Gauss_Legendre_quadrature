format longE
error(1) = abs(integral(@(x) x.^3+2.*x.^2, -1, 2) - integrate_using_Gauss_Legendre(-1, 2, 3, @(x) x.^3+2.*x.^2));
error(2) = abs(integral(@(x) sin(x)./x, 0, pi) - integrate_using_Gauss_Legendre(0, pi, 7, @(x) sin(x)./x));
error(3) = abs(integral(@(x) x.*exp(x), -20, -10.5) - integrate_using_Gauss_Legendre(-20, -10.5, 101, @(x) x.*exp(x)));