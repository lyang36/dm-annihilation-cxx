%%PLOTTING SCRIPT

%% Flux Profile
load output.txt
x=output(:,1);
y=output(:,2);

load out_core.txt
y1=out_core(:,2)

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
fdata1 = loglog(x*180/pi,y1, '-g', 'LineWidth',1.3);
fnfw = loglog(x*180/pi,vjnfw,'--r');
feinasto = loglog(x*180/pi,vjeinasto,'-.m');
fgnfw12 = loglog(x*180/pi,vjgnfw12);
fgnfw14 = loglog(x*180/pi,vjgnfw14);
ylabel('\Phi(GeV^2 cm^{-6} kpc sr^{-1})','FontSize',13);
xlabel('\theta/Deg');
title('Anglar Profile of Annahilation Flux','FontSize',13);
set(gca,'FontSize',12);
xlim([0.1 90])
legend('Regular Data', 'Cored Data', 'NFW', 'Einasto', 'GNFW1.2', 'GNFW1.4');


%% Density Profile
figure
load routput_host.txt
load rout_core1.txt
x=routput_host(:,1);
data=routput_host(:,2) / (routput_host(2,2));
dcore = rout_core1(:,2) / (rout_core1(200,2)) * data(200);
loglog(x*1000, data, '-k', 'LineWidth',1.3);
hold all
loglog(x*1000, dcore, '-.g', 'LineWidth',1.3);

rnfw = NNFW(x*1000,20);
%rnfw = rnfw / rnfw(2);
loglog(x*1000, rnfw, '--r');


reinasto = NEinasto(x*1000, 0.165, 20);
%reinasto = reinasto / reinasto(2);
loglog(x*1000, reinasto);

rgnfw = NGNFW(x*1000, 1.2, 20);
%rgnfw = rgnfw / rgnfw(2);
loglog(x*1000, rgnfw);

rpnfw = NPNFW(x*1000, 20);
%rgnfw = rgnfw / rgnfw(2);
loglog(x*1000, rpnfw);

ylabel('\rho (Normalized)','FontSize',13);
xlabel('r/kpc');
title('Density - Radial profile','FontSize',13);
set(gca,'FontSize',12);
xlim([x(2)*1000 100])
legend('Regular Data', 'CORED', 'NFW', 'Einasto', 'GNFW1.2', 'PNFW');



figure

loglog(x*1000, data.^2, '-k', 'LineWidth',1.3);
hold all
loglog(x*1000, dcore.^2, '-.g', 'LineWidth',1.3);

rnfw = NNFW(x*1000,20);
%rnfw = rnfw / rnfw(2);
loglog(x*1000, rnfw.^2, '--r');


reinasto = NEinasto(x*1000, 0.165, 20);
%reinasto = reinasto / reinasto(2);
loglog(x*1000, reinasto.^2);

rgnfw = NGNFW(x*1000, 1.2, 20);
%rgnfw = rgnfw / rgnfw(2);
loglog(x*1000, rgnfw.^2);

rpnfw = NPNFW(x*1000, 20);
%rgnfw = rgnfw / rgnfw(2);
loglog(x*1000, rpnfw.^2);

ylabel('\rho^2 (Normalized)','FontSize',13);
xlabel('r/kpc');
title('Density - Radial profile','FontSize',13);
set(gca,'FontSize',12);
xlim([x(2)*1000 100])
legend('Regular Data', 'CORED', 'NFW', 'Einasto', 'GNFW1.2', 'PNFW');


%% Density Profile
figure
load routall.txt
x=routall(:,1);
data=routall(:,2) / (routall(2,2));
loglog(x*1000, data, '-k', 'LineWidth',1.3);
hold all

rnfw = NNFW(x*1000,20);
%rnfw = rnfw / rnfw(2);
loglog(x*1000, rnfw, '--r');


reinasto = NEinasto(x*1000, 0.165, 20);
%reinasto = reinasto / reinasto(2);
loglog(x*1000, reinasto);

rgnfw = NGNFW(x*1000, 1.2, 20);
%rgnfw = rgnfw / rgnfw(2);
loglog(x*1000, rgnfw);

rpnfw = NPNFW(x*1000, 20);
%rgnfw = rgnfw / rgnfw(2);
loglog(x*1000, rpnfw);

ylabel('\rho (Normalized)','FontSize',13);
xlabel('r/kpc');
title('Density - Radial profile','FontSize',13);
set(gca,'FontSize',12);
xlim([x(2)*1000 100])
legend('Regular Data', 'NFW', 'Einasto', 'GNFW1.2', 'PNFW');


%%

i = 0;
filename = ['ann_nosop_sub', sprintf('%02i', i), '_21d.txt'];
d1 = textread(filename);
x = d1(:,1);
y1 = d1(:,2)/min(d1(:,2));
    
for i=1:size(x)
    vjnfw(i) = NJNFW(8, x(i), 20);
end
vjnfw = vjnfw / (min(vjnfw));

for i = 0:39
    figure('Visible','off')
    filename = ['ann_nosop_sub', sprintf('%02i', i), '_21d.txt'];
    d1 = textread(filename);
    x = d1(:,1);
    y1 = d1(:,2)/min(d1(:,2));
    
    filename = ['ann_sopv_sub', sprintf('%02i', i), '_21d.txt'];
    d2 = textread(filename);
    y2 = d2(:,2)/min(d2(:,2));
    
    filename = ['ann_sopv2_sub', sprintf('%02i', i), '_21d.txt'];
    d3 = textread(filename);
    y3 = d3(:,2)/min(d3(:,2));
    
    filename = ['decay_nosop_sub', sprintf('%02i', i), '_21d.txt'];
    d4 = textread(filename);
    y4 = d4(:,2)/min(d4(:,2));
    
    loglog(x * 180 / pi, y1, '-k');
    hold on
    
    loglog(x * 180 / pi, y2, '--r');

    loglog(x * 180 / pi, y3,'-.m');

    loglog(x * 180 / pi, y4, '-*b');  
    
    loglog(x * 180 / pi, vjnfw, '--g');  
    
    ylabel('\Psi (Normalized)');
    xlabel('\Theta (Deg)');
    title('Flux - Radial profile');
    xlim([0 20])
    ylim([min(y3), max(y3)])
    legend('No Correction', '1/v', '1/v^2', 'Decay', 'NFW-Annihilation');
    
    saveas(gcf,['angprof_', sprintf('%02i', i), '.png'],'png')
end




