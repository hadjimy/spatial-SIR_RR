function cubature_abscissas_plot(pts,cubature,xc,yc,delta)

n1 = pts;
n2 = 2*n1;

if strcmp(cubature,'elhay-kautsky')
    [~, x, y, ~, ~] = disk_rule(n1, n2, xc, yc, delta);
elseif strcmp(cubature,'gauss-legendre')
    [xp,~] = gauss_legendre_rule_compute(n1, 0, 1);
    [yp,~] = gauss_legendre_rule_compute(n2, 0, 1);
    [xx,yy] = meshgrid(xp,yp);
    x = delta.*xx.*cos(2*pi*yy);
    y = delta.*xx.*sin(2*pi*yy);
end

th = 0:pi/50:2*pi;
xunit = delta * cos(th) + xc;
yunit = delta * sin(th) + yc;

% plotting
fig = figure();
plot(xunit, yunit);
hold on
plot(x,y,'ko')
ticks = [-1 -delta/2 0 delta/2 1];
hc = get(fig,'children'); set(hc, 'fontsize', 14);
yticks(ticks); xticks(ticks);
axis('square')

filename = 'figures/'+string(cubature) + '_' + string(n1) + ...
    'x' + string(n2) + '.pdf';
saveas(gcf,filename)
