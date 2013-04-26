GNFW=@(x, Gamma, rs) 1./(x./rs).^Gamma./(1 + (x./rs)).^(3 - Gamma);
NFW=@(x, rs) 1./(x./rs)./(1 + (x./rs)).^2;
Einasto=@(x, Alpha, rs) exp(-2.0./Alpha.*((x./rs).^Alpha - 1));
PNFW= @(x, rs) (1./(x./rs).^1.5)./(1 + (x./rs)).^2;

NGNFW=@(x, Gamma, rs) GNFW(x, Gamma, rs)./GNFW(15, Gamma, 20).*0.01;
NNFW=@(x, rs) NFW(x, rs)./NFW(15, 20).*0.01;
NEinasto=@(x, Alpha, rs) Einasto(x, Alpha, rs)./Einasto(15, Alpha, 20).*0.01;
NPNFW = @(x, rs) PNFW(x, rs) ./ PNFW(15, 20).*0.01;

NJNFW = @(d, theta, rs) JNFW(d, theta, rs)./JNFW(d, 20*pi/180, rs)*2.91
NJEinasto = @(d, theta, Alpha, rs) JEinasto(d, theta, Alpha, rs) ...
    ./ JEinasto(d, 20*pi/180, Alpha, rs)*2.31
NJGNFW= @(d, theta, Alpha, rs) JGNFW(d, theta, Alpha, rs) ./ ...
    JGNFW(d, 20*pi/180, Alpha, rs)*2.91

NJPNFW = @(d, theta, rs) JPNFW(d, theta, rs)./JPNFW(d, 20*pi/180, rs)*2.91


NGNFW(8.5, 1.4, 20)
NEinasto(8.5, 0.165, 20)
NGNFW(8.5, 1.2, 20)
NGNFW(8.5, 1.4, 20)
NJNFW(8.50, 30/180*pi, 17)
NJEinasto(8.5, 30/180*pi, 0.165, 17)
NJGNFW(8.5, 30/180*pi, 1.2, 17)
NJGNFW(8.5, 30/180*pi, 1.4, 17)