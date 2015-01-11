function [ output ] = JNFW( d, theta, rs )
%%The annihilation of a NFW profile

NFW=@(x, rs) 1./(x./rs)./(1 + (x./rs)).^2;
NNFW=@(x, rs) NFW(x, rs)./NFW(8.5, 20).*0.22;

integrand = @(l) NNFW(sqrt((d - l .* cos(theta)).^2 + (l .* sin(theta)).^2), rs).^2;

output= 1/8.5 .* (1/0.3).^2 .* integral(integrand , 0, inf);

end

