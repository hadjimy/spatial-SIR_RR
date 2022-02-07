function fig1

quad_pts = [3,6,12];
cubature_type = {'elhay-kautsky','gauss-legendre'};
xc = 0;
yc = 0;
delta = 1;

for i = 1:length(cubature_type)
    cubature = cubature_type{i};
    
    for j = 1:length(quad_pts)
        pts = quad_pts(j);
        cubature_abscissas_plot(pts,cubature,xc,yc,delta)
    end
    
end
