function [order, error] = cubature_convergence(cubature,integrand)

quad_pts = [3,4,6,9,12];
colors = {[1 0 0],[0 0 1],[0 0 0],[0.75 0 0.75],[0 0.5 0]};

delta_size = 7; % length of delta vector
delta0 = 0.2; % The initial value of delta
delta = delta0./(2.^(0:delta_size-1)); % delta

% initialization
order = zeros(delta_size-1,length(quad_pts));
error = zeros(delta_size,length(quad_pts));
fig = figure();

for j = 1:length(quad_pts)
    for i = 1:delta_size
        aprx = cubature_int(delta(i),quad_pts(j),cubature,integrand);
        
        if strcmp(integrand,'initial')
            s = 0.1; % deviation for inital condition I0
            sol = 5000*(2*delta(i)-sqrt(2*pi)*s*erf(delta(i)/(sqrt(2)*s)));
        elseif  strcmp(integrand,'general')
            sol = (1/160).*pi.*(delta(i).*((-60)+delta(i).*...
                (30+delta(i).*((-20)+3.*(5+(-4).*delta(i)).*delta(i))))+...
                60.*log(1+delta(i)));
        end
        
        error(i,j) = abs(aprx-sol);
    end
    
    for i = 1:delta_size-1
        order(i,j) = ...
            real(log(error(i+1,j)/error(i,j))/log(delta(i+1)/delta(i)));
    end
    
    % plotting
    loglog(delta,error(:,j),'-d','color',colors{j},'LineWidth',1.2, ...
        'DisplayName',string(quad_pts(j))+' x '+string(2*quad_pts(j)))
    hold on
end

hc = get(fig,'children'); set(hc, 'fontsize', 12);
legend('-DynamicLegend','FontSize',14,'Location','northwest', ...
    'Interpreter','latex');
%set(hc,'fontsize',12,'XDir','reverse')
xlabel('$$\delta$$','FontSize',18,'Interpreter','latex');
ylabel('$$|error|$$','FontSize',18,'Interpreter','latex');
xticks(flip(delta));
%xlimits = xlim;
ylimits = ylim;
xlim([0.95*delta(end) 1.05*delta(1)])
ylim([ylimits(1)/10 10*ylimits(2)])

filename = append('figures/','convergence_',cubature,'_',integrand);
saveas(gcf,filename,'pdf')

end
