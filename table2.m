function table2

gridpts = 30;
quad_pts = 6; % 6 x 12 cubature points
tmax = 100;
discretization = 'elhay-kautsky';
%discretization = 'gauss-legendre';
interpol = 'linear';
method = 'FE';

delta = [0.025,0.05,0.075,0.1];
tau_pes = zeros(length(delta),1);
tau = zeros(length(delta),1);
tau_e = zeros(length(delta),1);

for i = 1:length(delta)
    fprintf('delta = %.3f\n',delta(i))
    [~,~,~,tau_pes(i),~] = spatial_sir_solve(gridpts,'pessimistic',tmax,...
        discretization,quad_pts,interpol,method,{'table_delta',delta(i)});
    [~,~,~,tau(i),~] = spatial_sir_solve(gridpts,'improved',tmax, ...
        discretization,quad_pts,interpol,method,{'table_delta',delta(i)});
    [~,~,~,tau_e(i),~] = spatial_sir_solve(gridpts,'adaptive',tmax, ...
        discretization,quad_pts,interpol,method,{'table_delta',delta(i)});
    fprintf('\n')
end

fprintf('delta \t\t tau_pes \t tau_pes/tau_e \t tau \t\t tau/tau_e \t tau_e \n')
for i = 1:length(delta)
    fprintf('%f\t %f\t %f\t %f\t %f\t %f\n', delta(i), tau_pes(i), ...
        tau_pes(i)/tau_e(i), tau(i), tau(i)/tau_e(i), tau_e(i))
end

fprintf('\nLatex output:\n')
for i = 1:length(delta)
    fprintf('$%.3f$ & $%.4f$ & $%.4f$ & $%.4f$ & $%.4f$ & $%.4f$ \\\\ \n',...
        delta(i), tau_pes(i), tau_pes(i)/tau_e(i), tau(i), ...
        tau(i)/tau_e(i), tau_e(i))
end

end
