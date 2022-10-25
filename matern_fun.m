function [y] = matern_fun(x)
    rho = 10;
    nu = 3.8;
    prefac = 2^(1-nu)/gamma(nu);
    xsc = sqrt(2*nu)/rho*x;
    y = prefac*xsc.^(nu).*besselk(nu,xsc);
end