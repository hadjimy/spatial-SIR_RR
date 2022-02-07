function fig2

cubature_types = {'elhay-kautsky','gauss-legendre'};
integrand = 'initial';

conv_res = containers.Map('KeyType','char','ValueType','any');

for k = 1:length(cubature_types)
    cubature = cubature_types{k};
    [order, error] = cubature_convergence(cubature,integrand);
    conv_res(cubature) = containers.Map({'order','error'},{order,error});
end

format long
for k = 1:length(cubature_types)
    results = conv_res(cubature_types{k});
    fprintf('%s:\n', cubature_types{k})
    fprintf('  order \n')
    disp(results('order'))
    fprintf('  error \n')
    disp(results('error'))
end

end
