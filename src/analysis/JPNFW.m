function [ output ] = JPNFW( d, theta, rs )
%%The annihilation of a NFW profile

    PNFW= @(x, rs) 1./(x./rs).^{1.5}/(1 + (x./rs)).^2;
    integrand = @(l) PNFW(sqrt((d - l .* cos(theta)).^2 + (l .* sin(theta)).^2), rs).^2;

    output= 1/8.5 .* (1/0.3).^2 .* integral(integrand , 0, inf);

end

