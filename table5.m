function table5
% Solve the SIR problem and save solution

gridpts = 30;
quad_pts = 6; % 6 x 12 cubature points
tbound = 'improved';
%discretizations = {'elhay-kautsky','gauss-legendre'};
discretizations = {'gauss-legendre'};
interpolations = {'linear','makima'};
methods = {'FE','SSPRK33','SSPRK104'};

tAll = zeros(length(methods),length(discretizations),length(interpolations));
tSetup = zeros(size(tAll));
tSpatial = zeros(size(tAll));
tMethod = zeros(size(tAll));

t = 50;
for i = 1:length(methods)
    for k = 1:length(interpolations)
        for j = 1:length(discretizations)
            tStart = tic;
            [~,~,~,~,~,CPUtimes] = spatial_sir_solve(gridpts,tbound,t, ...
                discretizations{j},quad_pts,interpolations{k}, ...
                methods{i},{'default'},NaN,'checks_off');
            tAll(i,j,k) = toc(tStart);
            tSetup(i,j,k) = CPUtimes(1);
            tSpatial(i,j,k) = CPUtimes(2);
            tMethod(i,j,k) = CPUtimes(3);
        end
    end
end

%% Output results
for i = 1:length(methods)
    fprintf('Method: %s \n',methods{i})
    fprintf(['\t Interpolation \t Quadrature \t Time-step setup' ...
        '\t Spatial discretization \t Method \t Total\n'])
    for k = 1:length(interpolations)
        for j = 1:length(discretizations)
            fprintf('\t %s \t %s\t %f\t\t %f\t\t\t %f\t %f\n', ...
                interpolations{k}, discretizations{j}, tSetup(i,j,k), ...
                tSpatial(i,j,k), tMethod(i,j,k), tAll(i,j,k))
        end
    end
end

fprintf('\n\n')

for i = 1:length(methods)
    fprintf('Method: %s \n',methods{i})
    fprintf(['\t Interpolation \t Quadrature \t Spatial discretization '...
        '\t Total \n'])
    for k = 1:length(interpolations)
        for j = 1:length(discretizations)
            fprintf('\t %s \t %s\t %f\t\t\t %f\n', ...
                interpolations{k}, discretizations{j}, ...
                tSetup(i,j,k)+tSpatial(i,j,k), tAll(i,j,k))
        end
    end
end

fprintf('\nLatex output:\n')
interpolations = {'bilinear','makima'};
for i = 1:length(methods)
    for k = 1:length(interpolations)
        for j = 1:length(discretizations)
            if k == 1
                fprintf('\\multirow{2}{*}{%-8s%s', methods{i})
                fprintf('} & ')
            else
                fprintf('\t\t\t  & ')
            end
            fprintf('%s \t& $%6.3f$ & $%6.3f$ \\\\ \n', ...
                interpolations{k}, tSetup(i,j,k)+tSpatial(i,j,k), tAll(i,j,k))
        end
    end
end

end
