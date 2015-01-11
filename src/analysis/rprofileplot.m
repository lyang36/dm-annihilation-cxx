function [x y] = rprofileplot(filename, outputfile, setxlim = false, setylim=false, sxlim=[-1, 1], sylim=[-1,1])
%% Plot R-profile
fid = fopen(filename);
fgetl(fid);
fgetl(fid);
fgetl(fid);
A = fscanf(fid, '%f %f', [2, inf]);
x = (A(1,:))' / pi * 180;
y = (A(2,:))';

y = y / max(y);

figure('Visible','off');
%figure;
loglog(x, y, 'k');
ylabel('\psi (Normalized)','FontSize',13);
xlabel('angular radius (degree)');
title('Angular profile','FontSize',13);
set(gca,'FontSize',12);
if (setxlim)
    xlim(sxlim)
else
    xlim([min(x), max(x)]);
end
if (setylim)
    ylim(sylim)
else
    ylim([min(y), max(y)]);
end
saveas(gcf,[outputfile, '.png'],'png');

end


