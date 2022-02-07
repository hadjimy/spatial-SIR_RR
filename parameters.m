function [a,b,c,alpha,beta,delta,tol] = parameters(gridpts,output)

tol = 10^(-12);

if strcmp(output{1},'default') || strcmp(output{1},'plot_fig')
    a = 100; b = 0.05; c = 0.01; delta = 0.05;
elseif strcmp(output{1},'plot_error')
    a = 200; b = 0.01; c = 0.01; delta = 0.1;
elseif strcmp(output{1},'table_a')
    a = output{2}; b = 0.05; c = 0.01; delta = 0.05;
elseif strcmp(output{1},'table_delta')
    a = 100; b = 0.05; c = 0.01; delta = output{2};
end

% load air currents over North America from MATLAB's dataset
load('wind','u','v')
u = [u(:,:,1), u(:,:,2); u(:,:,3), u(:,:,4)];
v = [v(:,:,1), v(:,:,2); v(:,:,3), v(:,:,4)];

% find total points of grid
h = 1/(gridpts-1);
gh = ceil(delta/h);
N = gridpts + 2*gh;
[X,Y] = meshgrid(-gh*h:h:1+gh*h,1+gh*h:-h:-gh*h);

% directional components of wind vector at each grid point  
U = u(1:N,1:N,1);
V = v(1:N,1:N,1);
Z = sqrt(U.^2 + V.^2);

% reduce magnitude
for i = 1:N
    for j = 1:N
        U(i,j) = U(i,j)/norm([U(i,j), V(i,j)]);
        V(i,j) = V(i,j)/norm([U(i,j), V(i,j)]);
    end
end

% % normalize U and Y
% U = zeros(N);
% V = zeros(N);
% for i = 1:N
%     for j = 1:N
%         U(i,j) = u(i,j)/norm([u(i,j), v(i,j)]);
%         V(i,j) = v(i,j)/norm([u(i,j), v(i,j)]);
%     end
% end

if strcmp(output{1},'plot_fig')
    contourf(X,Y,Z,10,'LineColor','none');
    hold on
    fig = quiver(X,Y,U,V,'Color','w');
    axis([0 1 0 1])
    hc = get(fig,'children');
    set(hc,'FontSize',14)
    filename = append('figures/wind_profile');
    saveas(fig,filename,'pdf')
end

% find angles and magnitude of wind at each grid point 
alpha = zeros(N,N);
beta = zeros(N,N);
for i = 1:N
    for j = 1:N
        if U(i,j) ~= 0
            angle = atan(abs(V(i,j)/U(i,j)));
        else
            if V(i,j) > 0
                alpha(i,j) = pi/2;
            elseif V(i,j) < 0
                alpha(i,j) = 3*pi/2;
            else
                alpha(i,j) = 0;
            end
        end
        
        if U(i,j) > 0 && V(i,j) >= 0
            alpha(i,j) = angle;
        elseif U(i,j) < 0 && V(i,j) >= 0
            alpha(i,j) = pi - angle;
        elseif U(i,j) < 0 && V(i,j) <= 0
            alpha(i,j) = pi + angle;
        elseif U(i,j) > 0 && V(i,j) <= 0
            alpha(i,j) = 2*pi - angle;
        end
        
        beta(i,j) = norm([U(i,j), V(i,j)]);
    end
end

end
