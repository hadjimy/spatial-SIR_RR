x = -2:0.2:2;
y = x';
z = x .* exp(-x.^2 - y.^2);
%z = peaks(x,y);
[px,py] = gradient(z);

[X,Y] = meshgrid(x,y);
Z = X.*exp(-X.^2 - Y.^2);
[U,V,W] = surfnorm(X,Y,z);

figure(1)
surfc(x,y,z);
hold on
zz = min(z,[],'All')*ones(size(z));
quiver3(X,Y,zz,U,V,zeros(size(Z)));


figure(2)
sc = contour(x,y,z);
hold on
qv = quiver(x,y,-px,-py);



