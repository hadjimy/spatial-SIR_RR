function Q = cubature_int(delta,quad_pts,cubature,integrand_type)
% Calculates the weighted integral of the function 'integrand' 
% (defined at the end of this code) at point (0,0) over a circle with
% radius 'delta'.

% problem parameters
a = 100; beta = 1; alpha = 0;

n1 = quad_pts;
n2 = 2*n1;
n = n1*n2;

if strcmp(cubature,'elhay-kautsky')
    [w1, ~, ~, r, theta] = disk_rule(n1, n2, 0, 0, delta);
    
    quad_pts = {r,theta};
    w = {w1,[]};
    
    integrand=zeros(n,1);
    for j = 1:n
        integrand(j)=fun(r(j),theta(j),integrand_type,a,alpha,beta,delta);
    end
    
elseif strcmp(cubature,'gauss-legendre')
    %[z, w1] = quadrature_pts(n1);
    
    %[z, w1] = line_ncc_rule ( n1, 0, 1 ); %% negative weights
    %[z, w1] = line_nco_rule ( n1, 0, 1 ); %% negative weights
    %[z, w1] = ccn_rule ( n1, 0, 1, 'test' ); %% negative weights
    %[z, w1] = clenshaw_curtis_rule ( n1, 0, (n1-1)/n1, 'test' );
    
    %[z, w1] = patterson_rule ( n1, 0, 1, 'test' ); %% same as legendre
    %[z, w1] = gegenbauer_rule ( n1, 0, 0, 1, 'test' ); %% same as legendre
    %[z, w1] = jacobi_rule ( n1, 0, 0, 0, 1, 'test' ); %% same as legendre
    %[z, w1] = legendre_rule ( n1, 0, 1, 'test' );
    [z, w1] = gauss_legendre_rule_compute (n1, 0, 1 );
    
    %[eta, w2] = quadrature_pts(n2);
    
    %[eta, w2] = line_ncc_rule ( n2,  0, (n2-1)/n2 ); %% negative weights
    %[eta, w2] = line_nco_rule ( n2, 0, 1 ); %% negative weights
    %[eta, w2] = ccn_rule ( n2, 0, 1, 'test' ); %% negative weights
    %[eta, w2] = clenshaw_curtis_rule ( n2, 0, (n2-1)/n2, 'test' );
    
    %[eta, w2] = patterson_rule ( n2, 0, 1, 'test' ); %% same as legendre
    %[eta, w2] = gegenbauer_rule ( n2, 0, 0, 1, 'test' ); %% same as legendre
    %[eta, w2] = jacobi_rule ( n2, 0, 0, 0, 1, 'test' ); %% same as legendre
    %[eta, w2] = legendre_rule ( n2, 0, 1, 'test' );
    [eta, w2] = gauss_legendre_rule_compute (n2, 0, 1 );
    
    quad_pts = {z, eta};
    w = {w1,w2};
    
    % The values which will be written inside the sum of the numerical integral:
    integrand = zeros(n1,n2);
    
    for j = 1:n1
        for k = 1:n2
            r = z(j)*delta;
            theta = 2*pi*eta(k);
            integrand(j,k)=fun(r,theta,integrand_type,a,alpha,beta,delta);
        end
    end
    
end

Q = numint(delta,integrand,quad_pts,w,cubature);

end

function sum = numint(delta, integrand, pts, w, cubature)

sum = 0;
if strcmp(cubature,'elhay-kautsky')
    w = w{1};
    
    for j=1:length(w)
        sum = sum + w(j)*integrand(j);
    end
    sum = pi*delta^2*sum;
elseif strcmp(cubature,'gauss-legendre')
    z = pts{1};
    eta = pts{2};
    w_z = w{1};
    w_eta = w{2};
    
    for j = 1:length(z)
        for k = 1:length(eta)
           sum = sum + w_z(j)*w_eta(k)*2*pi*delta^2*z(j)*integrand(j,k);
        end
    end
end

end

function val = f1(r, a, delta)
% Determines the values of function f1.
val = a*(-r + delta);
end

function val = f2(theta, alpha, beta)
% Determines the values of function f2.
val = beta*sin(theta + alpha) + beta;
end

function val = fun(r, theta, integrand_type, a, alpha, beta, delta)

x = r*cos(theta);
y = r*sin(theta);

if strcmp(integrand_type,'initial')
    s=1/10; % deviation
    %s=delta;
    I = 100/(2*pi*s^2)*exp(-((x).^2/s^2+(y).^2/s^2)/2);
    val = f1(r,a,delta)*f2(theta,alpha,beta)*I;

elseif strcmp(integrand_type,'general')
    val = 5.*(1+r).^(-1).*(5.*r.*cos(theta)+(-2/3).*r.^3.*cos(theta).^3+...
       (-1/10).*r.^4.*sin(theta).^4+2.*r.^7.*cos(theta).^2.*sin(theta).^5);
end

end
