function fig4

Lp = 2;
gridpts = 30;
tmax = 50;
method = 'SSPRK104';

interpol = {'linear','makima','spline'};
discretization = {'elhay-kautsky','gauss-legendre'};
discretization_legend = {'Elhay-Kautsky','Gauss-Legendre'};
quad_pts = [linspace(4,10,7),20];
nodes = string(quad_pts(1:end-1))+' x '+string(2*quad_pts(1:end-1));
marks = append([":"; "-"],["ko" "bd" "rs"]);
error = containers.Map('KeyType','char','ValueType','any');
fig = figure();

h = 1/(gridpts-1);
for i = 1:length(discretization)
    for j = 1:length(interpol)
        fprintf('discretization: %s\n', discretization{i})
        fprintf('interpolation: %s\n', interpol{j})
        err = zeros(length(quad_pts)-1,1);
        
        [S,I,R,~,~] = spatial_sir_solve(gridpts,'improved',tmax, ...
            discretization{i},quad_pts(end),interpol{j},method, ...
            {'default'});
        [m,n] = size(S);
        ref_sol = [reshape(S,m*n,[]);reshape(I,m*n,[]);reshape(R,m*n,[])];
        
        for k = 1:length(quad_pts)-1
            [S,I,R,~,~] = spatial_sir_solve(gridpts,'improved',tmax, ...
                discretization{i},quad_pts(k),interpol{j},method, ...
                {'default'});
            [m,n] = size(S);
            sol = [reshape(S,m*n,[]);reshape(I,m*n,[]);reshape(R,m*n,[])];
            err(k) = (h^2)^(1/Lp)*norm(sol-ref_sol,Lp);
        end
        fprintf('\n')
        
        if strcmp(interpol{j},'linear')
            interpol_name = 'bilinear';
        else
            interpol_name = interpol{j};
        end
        key = append(discretization_legend{i},'; ',interpol_name);
        error(key) = err;
        
        % plotting
        loglog(quad_pts(1:end-1),error(key),marks{i,j},'LineWidth',1.2, ...
            'MarkerSize',8,'DisplayName',key)
        hold on
        
    end
end

hc = get(fig,'children'); set(hc, 'fontsize', 12);
legend('-DynamicLegend','FontSize',16,'Location','SouthWest', ...
    'Interpreter','latex');
xticks(quad_pts(1:end-1))
xticklabels(nodes)
ylimits = ylim;
ylim([ylimits(1)/100 ylimits(2)])
xlabel('$$cubature \, nodes$$','FontSize',18,'Interpreter','latex');
ylabel('$$|error|$$','FontSize',18,'Interpreter','latex');

filename = append('figures/','cubature_convergence');
saveas(gcf,filename,'pdf')

end
