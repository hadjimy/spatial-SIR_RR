function table1

gridpts = 30;
quad_pts = 6; % 6 x 12 cubature points
tmax = 100;
discretization = 'elhay-kautsky';
%discretization = 'gauss-legendre';
interpol = 'linear';
method = 'FE';

a = [50,100,250,500];
tau_pes = zeros(length(a),1);
tau = zeros(length(a),1);
tau_e = zeros(length(a),1);

for i = 1:length(a)
    fprintf('a = %d\n',a(i))
    [~,~,~,tau_pes(i),~] = spatial_sir_solve(gridpts,'pessimistic',tmax,...
        discretization,quad_pts,interpol,method,{'table_a',a(i)});
    [~,~,~,tau(i),~] = spatial_sir_solve(gridpts,'improved',tmax, ...
        discretization,quad_pts,interpol,method,{'table_a',a(i)});
    [~,~,~,tau_e(i),~] = spatial_sir_solve(gridpts,'adaptive',tmax, ...
        discretization,quad_pts,interpol,method,{'table_a',a(i)});
    fprintf('\n')
end

fprintf('a \t\t tau_pes \t tau_pes/tau_e \t tau \t\t tau/tau_e \t tau_e \n')
for i = 1:length(a)
    fprintf('%f\t %f\t %f\t %f\t %f\t %f\n', a(i), tau_pes(i), ...
        tau_pes(i)/tau_e(i), tau(i), tau(i)/tau_e(i), tau_e(i))
end

fprintf('\nLatex output:\n')
for i = 1:length(a)
    fprintf('$%d$ \t & $%.4f$ & $%.4f$ & $%.4f$ & $%.4f$ & $%.4f$ \\\\ \n',...
        a(i), tau_pes(i), tau_pes(i)/tau_e(i), tau(i), tau(i)/tau_e(i), ...
        tau_e(i))
end

end
