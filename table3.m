function table3

Lp = 2;
gridpts = 30;
quad_pts = 6; % 6 x 12 cubature points
tmax = 50;
discretization = 'gauss-legendre';
interpol = 'spline';
methods = {'num_integral','FE'};
error = containers.Map('KeyType','char','ValueType','any');
order = containers.Map('KeyType','char','ValueType','any');

%% Convergence test
fprintf(['Running code with %s method and %s time-step restriction to ' ...
    'find largest allowed time step ...\n'],'FE', 'adaptive')

[~,~,~,largest_timestep,~] = spatial_sir_solve(gridpts,'adaptive',tmax, ...
        discretization,quad_pts,interpol,'FE',{'default'});
%largest_timestep = floor(largest_timestep);
largest_timestep = num2str(largest_timestep,'%.2f');
largest_timestep = str2double(largest_timestep(1:end-1));
timesteps = largest_timestep*2.^(0:-1:-8);

fprintf(['Time steps to use: ', repmat('%g, ', 1, numel(timesteps)-1), ...
    '%g\n\n'],timesteps)
fprintf('Computing solution with the following methods:\n')

h = 1/(gridpts-1);
for i = 1:length(methods)
    
    if strcmp(methods{i},'num_integral')
        discretization = 'num_integral';
    else
        discretization = 'gauss-legendre';
    end
    
    disp(methods{i})
    [Sref,Iref,Rref,~,~] = spatial_sir_solve(gridpts,'given',tmax, ...
        discretization,quad_pts,interpol,methods{i},{'default'}, ...
        timesteps(end)/2);
    [m,n] = size(Sref);
    ref_sol = [reshape(Sref,m*n,[]);reshape(Iref,m*n,[]); ...
        reshape(Rref,m*n,[])];
    
    err = zeros(length(timesteps),1);
    ord = zeros(length(err)-1,1);
    
    for j = 1:length(timesteps)
        [S,I,R,~,~] = spatial_sir_solve(gridpts,'given',tmax, ...
        discretization,quad_pts,interpol,methods{i},{'default'}, ...
        timesteps(j));
        [m,n] = size(S);
        sol = [reshape(S,m*n,[]);reshape(I,m*n,[]);reshape(R,m*n,[])];
        
        err(j) = (h^2)^(1/Lp)*norm(sol-ref_sol,Lp);
    end
    
    for k = 1:length(ord)
        ord(k) = log(err(k+1)/err(k))/ ...
            log(timesteps(k+1)/timesteps(k));
    end
    
    error(methods{i}) = err;
    order(methods{i}) = ord;
    
end

for i = 1:length(methods)
    mthd_error =  error(methods{i});
    mthd_order =  order(methods{i});
    fprintf('Method: %s \n',methods{i})
    fprintf('timestep \t error \t\t order \n')
    fprintf('%f\t %.2e\n', timesteps(1), mthd_error(1))
    for j = 2:length(timesteps)
        fprintf('%f\t %.2e\t %.2f\n', timesteps(j), mthd_error(j),...
        mthd_order(j-1))
    end
end

fprintf('\nLatex output:\n')
for i = 1:length(methods)
    mthd_error =  error(methods{i});
    mthd_order =  order(methods{i});
    
    str_error = sprintf('%.2e',mthd_error(1));
    fprintf('Method: %s \n',methods{i})
    fprintf('$%.4f$ & $%s \\times 10^{%s}$ & \\\\ \n', timesteps(1), ...
        str_error(1:4), str_error(end-2:end))
    for j = 2:length(timesteps)
        str_error = sprintf('%.2e',mthd_error(j));
        fprintf('$%.4f$ & $%s \\times 10^{%s}$ & $%.2f$ \\\\ \n', ...
            timesteps(j), str_error(1:4), str_error(end-2:end), ...
            mthd_order(j-1))
    end
    fprintf('\n')
end

end
