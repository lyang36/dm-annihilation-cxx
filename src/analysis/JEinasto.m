function [ output ] = JEinasto(d, theta, Alpha, rs)
%%The annihilation of a NFW profile

    Einasto=@(x, Alpha, rs) exp(-2.0./Alpha.*((x./rs).^Alpha - 1));
    NEinasto=@(x, Alpha, rs) Einasto(x, Alpha, rs)./Einasto(8.5, Alpha, 20).*0.22;

    integrand = @(l)  NEinasto(sqrt((d - l .* cos(theta)).^2 + (l .* sin(theta)).^2), Alpha, rs).^2;

    output= 1/8.5 .* (1/0.3).^2 .* integral(integrand , 0, inf);

end

