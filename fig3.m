function fig3
% Solve the SIR problem and save solution

gridpts = 30;
quad_pts = 6; % 6 x 12 cubature points
tbound = 'improved';
discretization = 'gauss-legendre';
interpol = 'makima';
method = 'SSPRK104';

t = 50;
[~,~,~,tau,fig] = spatial_sir_solve(gridpts,tbound,t,discretization, ...
    quad_pts,interpol,method,{'plot_fig'});
hc = get(fig,'children');
hc(2).Title.String = ' '; % remove title
set(hc,'FontSize',8)
filename = append('figures/','sol_',discretization,'_t=',string(t));
saveas(fig,filename,'pdf')

close all
t = 1000;
[~,~,~,~,fig] = spatial_sir_solve(gridpts,tbound,t,discretization, ...
    quad_pts,interpol,method,{'plot_fig'});
hc = get(fig,'children');
hc(2).Title.String = ' '; % remove title
set(hc,'FontSize',10)
filename = append('figures/','sol_',discretization,'_t=',string(t));
saveas(fig,filename,'pdf')

fprintf('time step = %f\n',tau)
end
