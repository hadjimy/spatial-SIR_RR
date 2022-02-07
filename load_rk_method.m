function [r,delta,beta,K,c] = load_rk_method(method)
% Load Shu-Osher and Butcher coefficients of Runge-Kutta methods
%
% Input variables:
%   method      --- explicit methods:
%                       'FE','MTE22','RK33','RK44',,'SSPRKsp',
%                       where s = stages, p = order
%               --- implicit methods:
%                       'BE','TM','IRK22u','IRK33u'
%
% Output variables:
%	r           --- SSP coefficient
%   delta,beta  --- Shu-Osher coefficients
%   K, c        --- Butcher coefficients
%==========================================================================

switch method
    case 'FE'
        % delta
        delta = sparse(2,1);
        delta(1,1) = 1.000000000000000;
        % beta
        beta = sparse(2,2);
        beta(2,1) = 1.000000000000000;
        %r
        r = 1.000000000000000;
    case 'MTE22'
        % delta
        delta = sparse(3,1);
        delta(1,1) = 1.000000000000000;
        delta(2,1) = 0.666666666666667;
        delta(3,1) = 0.625000000000000;
        % beta
        beta = sparse(3,3);
        beta(2,1) = 0.666666666666667;
        beta(3,2) = 0.750000000000000;
        %r
        r = 0.500000000000000;
    case 'RK33'
        % delta
        delta = sparse(4,1);
        delta(1,1) = 1.000000000000000;
        delta(2,1) = 1.000000000000000;
        delta(3,1) = 1.000000000000000;
        delta(4,1) = 1.000000000000000;
        % beta
        beta = sparse(4,4);
        beta(2,1) = 0.500000000000000;
        beta(4,1) = 0.166666666666667;
        beta(3,2) = 2.000000000000000;
        beta(4,2) = 0.666666666666667;
        beta(4,3) = 0.166666666666667;
        %r
        r = 0.000000000000000;
    case 'RK44'
        % delta
        delta = sparse(5,1);
        delta(1,1) = 1.000000000000000;
        delta(2,1) = 1.000000000000000;
        delta(3,1) = 1.000000000000000;
        delta(4,1) = 1.000000000000000;
        delta(5,1) = 1.000000000000000;
        % beta
        beta = sparse(5,5);
        beta(2,1) = 0.500000000000000;
        beta(5,1) = 0.166666666666667;
        beta(3,2) = 0.500000000000000;
        beta(5,2) = 0.333333333333333;
        beta(4,3) = 1.000000000000000;
        beta(5,3) = 0.333333333333333;
        beta(5,4) = 0.166666666666667;
        %r
        r = 0.000000000000000;
    case 'SSPRK22'
        % delta
        delta = sparse(3,1);
        delta(1,1) = 1.000000000000000;
        delta(3,1) = 0.500000000000000;
        % beta
        beta = sparse(3,3);
        beta(2,1) = 1.000000000000000;
        beta(3,2) = 0.500000000000000;
        %r
        r = 1.000000000000000;
    case 'SSPRK33'
        % delta
        delta = sparse(4,1);
        delta(1,1) = 1.000000000000000;
        delta(3,1) = 0.750000000000000;
        delta(4,1) = 0.333333333333333;
        % beta
        beta = sparse(4,4);
        beta(2,1) = 1.000000000000000;
        beta(3,2) = 0.250000000000000;
        beta(4,3) = 0.666666666666667;
        %r
        r = 1.000000000000000;
    case 'SSPRK104'
        % delta
        delta = sparse(11,1);
        delta(1,1) = 1.000000000000000;
        delta(4,1) = 0.000000000000008;
        delta(6,1) = 0.599999999999999;
        delta(7,1) = 0.000000000000001;
        delta(10,1) = 0.000000000000000;
        delta(11,1) = 0.039999999999999;
        % beta
        beta = sparse(11,11);
        beta(2,1) = 0.166666666666667;
        beta(3,2) = 0.166666666666667;
        beta(4,3) = 0.166666666666666;
        beta(5,4) = 0.166666666666667;
        beta(6,5) = 0.066666666666668;
        beta(10,5) = 0.000000000000003;
        beta(11,5) = 0.060000000000000;
        beta(7,6) = 0.166666666666667;
        beta(8,7) = 0.166666666666667;
        beta(9,8) = 0.166666666666667;
        beta(11,8) = 0.000000000000001;
        beta(10,9) = 0.166666666666664;
        beta(11,10) = 0.100000000000000;
        %r
        r = 5.999999999999984;
    case 'BE'
        % delta
        delta = sparse(2,1);
        delta(1,1) = 0.090909090909091;
        delta(2,1) = 0.090909090909092;
        % beta
        beta = sparse(2,2);
        beta(1,1) = 0.090909090909091;
        beta(2,1) = 0.090909090909092;
        %r
        r = 10.000000000000000;
    case 'TM'
        % delta
        delta = sparse(3,1);
        delta(1,1) = 1.000000000000000;
        % beta
        beta = sparse(3,3);
        beta(2,1) = 0.250000000000000;
        beta(3,1) = 0.250000000000000;
        beta(2,2) = 0.250000000000000;
        beta(3,2) = 0.250000000000000;
        %r
        r = 2.000000000000000;
end

I = speye(size(beta,1));
alpha = r*beta;
K = (I-alpha)\beta;
% abcissa
c = sum(K(1:end-1,1:end-1),2);

end
