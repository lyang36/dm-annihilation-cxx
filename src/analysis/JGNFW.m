function [ output ] = JGNFW( d, theta, Gamma, rs)
%%The annihilation of a NFW profile

    GNFW=@(x, Gamma, rs) 1./(x./rs).^Gamma./(1 + (x./rs)).^(3 - Gamma);

    NGNFW=@(x, Gamma, rs) GNFW(x, Gamma, rs)./GNFW(8.5, Gamma, 20).*0.22;

    integrand = @(l)  NGNFW(sqrt((d - l .* cos(theta)).^2 + (l .* sin(theta)).^2), Gamma, rs).^2;

    output= 1/8.5 .* (1/0.3).^2 .* integral(integrand , 0, inf);

end

