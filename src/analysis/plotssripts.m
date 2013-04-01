%%PLOTTING SCRIPT

load output.txt
x=output(:,1);
y=output(:,2);


vjnfw = x;
vjeinasto = x;
vjgnfw12 = x;
vjgnfw14 = x;
for i=1:size(x)
    vjnfw(i) = NJNFW(8, x(i), 20);
    vjeinasto(i) = NJEinasto(8, x(i), 0.165, 20);
    vjgnfw12(i) = NJGNFW(8, x(i), 1.2, 20);
    vjgnfw14(i) = NJGNFW(8, x(i), 1.4, 20);
end

figure
fdata = loglog(x*180/pi,y, '-k', 'LineWidth',1.3);
hold all
fnfw = loglog(x*180/pi,vjnfw,'--r');
feinasto = loglog(x*180/pi,vjeinasto,'-.m');
fgnfw12 = loglog(x*180/pi,vjgnfw12);
fgnfw14 = loglog(x*180/pi,vjgnfw14);
ylabel('\Phi(GeV^2 cm^{-6} kpc sr^{-1})','FontSize',13);
xlabel('\theta/Deg');
title('Anglar Profile of Annahilation Flux','FontSize',13);
set(gca,'FontSize',12);
xlim([0.1 90])
legend('Regular Data','NFW', 'Einasto', 'GNFW1.2', 'GNFW1.4');