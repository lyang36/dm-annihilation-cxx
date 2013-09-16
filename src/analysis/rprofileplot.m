function [x y] = rprofileplot(filename, outputfile)
%% Plot R-profile
fid = fopen(filename);
fgetl(fid);
fgetl(fid);
fgetl(fid);
A = fscanf(fid, '%f %f', [2, inf]);
x = (A(1,:))'
y = (A(2,:))'

figure('Visible','off');
%figure;
loglog(x/max(x), y, 'k');
ylabel('\rho','FontSize',13);
xlabel('r/r_0');
title('Density - Radial profile','FontSize',13);
set(gca,'FontSize',12);
xlim([0, 1]);
ylim([min(y), max(y)]);
saveas(gcf,[outputfile, '.png'],'png')

end


