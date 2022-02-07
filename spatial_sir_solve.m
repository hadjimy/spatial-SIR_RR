function [S,I,R,timestep,fig,CPUtimes] = spatial_sir_solve(gridpts,tbound,...
    tmax,discretization,quad_pts,interpol,method,output,tau_fixed,checks)

% Choice of timestep bound:
% tbound = 'pessimistic'
%        = 'improved'
%        = 'adaptive'
%        = 'given'
%
% Choice of discretization (using analytical solution or cubature rule):
% discretization = 'num_integral' (approximation of analytical solution)
%                = '4point' (the 4 point method)
%                = 'elhay-kautsky'
%                = 'gauss-legendre'
%
% Choice of interpolation:
% interpol = 'linear'
%          = 'cubic'
%          = 'spline'
%          = 'makima'
%
% Choice of output:
% output = 'plot_fig'
%        = 'plot_error'
%        = 'table_a'
%        = 'table_delta'
%
% Choice of quadrature points (quad_pts):
%   N_r = [3,4,6,9,13];
%   N_t = 2*N_r
%
% Example:
% [S,I,R,new_tau] = spatial_sir_solve(30,'improved',80, ...
%   'elhay-kautsky',13,'linear','SSPRK104','plot_fig');

%% Setup
if nargin > 9
    property_checks = checks;
else
    property_checks = 'checks_on';
end
% Load (fixed) parameters
[a,b,c,alpha,beta,delta,tol] = parameters(gridpts,output);

% Spatial grid
h = 1/(gridpts-1);
gh = ceil(delta/h);
N = gridpts + 2*gh;
[X,Y] = meshgrid(0:h:1,1:-h:0);
[XX,YY] = meshgrid(-gh*h:h:1+gh*h,1+gh*h:-h:-gh*h);

% Initial conditions
% prealocation
S = zeros(N); 
I = zeros(N);
R = zeros(N);
Snew = zeros(gridpts);
Inew = zeros(gridpts);
Rnew = zeros(gridpts);

S0 = 20; % initial susceptibles
s = 0.1; % deviation
I0 = 1/(2*pi*s^2);
%S(1+gh:N-gh,1+gh:N-gh) = S0*ones(gridpts);
I(1+gh:N-gh,1+gh:N-gh) = I0*exp(-((X-1/2).^2/s^2+(Y-1/2).^2/s^2)/2);
                         % + I0*exp(-((X-1/4).^2/s^2+(Y-3/4).^2/s^2)/4);
S(1+gh:N-gh,1+gh:N-gh) = I0*ones(gridpts) - I(1+gh:N-gh,1+gh:N-gh);

initsum = S(1+gh:N-gh,1+gh:N-gh) + I(1+gh:N-gh,1+gh:N-gh) + ...
    R(1+gh:N-gh,1+gh:N-gh);
%M = max(max(initsum));
M = max(initsum,[],'all');

%f1 = @(r) a*(-r+delta);
%f2 = @(t) beta*sin(t+alpha)+beta;
k1 = a*(2-sqrt(2))*delta/2; % maximum of f1;
k2 = 2*max(beta,[],'all'); % maximum of f2;

%% Load time-stepping method
if ~strcmp(discretization,'num_integral')
    [sspcoef,RK_delta,RK_beta] = load_rk_method(method);
    if sspcoef == 0.
        sspcoef = 1; % avoid a zero time step if RK method has sspcoef=0.
    end
    RK_delta = full(RK_delta);
    RK_beta = full(RK_beta);
    s = length(RK_delta)-1;
end


if strcmp(tbound,'adaptive') &&  ~strcmp(method,'FE')
    disp('The adaptive timestep only works for Forward Euler!')
    return
end

%% Calculate Tbar
tStart = tic; % start clock 
Tbar_all = zeros(gridpts);
for k = 1+gh:N-gh
    for l = 1+gh:N-gh
        if strcmp(discretization,'4point')
            n = 4;
            maxw = pi*delta^2/4;
            Tbar_all = maxw*k1*M*(f2(pi/4,alpha(k,l),beta(k,l)) + ...
                f2(3*pi/4,alpha(k,l),beta(k,l)) + ...
                f2(5*pi/4,alpha(k,l),beta(k,l)) + ...
                f2(7*pi/4,alpha(k,l),beta(k,l)));
        else
            n1 = quad_pts;
            n2 = 2*n1;
            n = n1*n2;
            % get quadrature points
            if strcmp(discretization,'elhay-kautsky')
                [w, ~, ~, r, theta] = disk_rule(n1, n2, 0, 0, delta);
            elseif strcmp(discretization,'gauss-legendre')
                [z, w1] = gauss_legendre_rule_compute(n1, 0, 1);
                [eta, w2] = gauss_legendre_rule_compute(n2, 0, 1);
            end
            
            %The construction:
            sum = 0;
            maxw = 0;
            if strcmp(discretization,'elhay-kautsky')
                for i=1:n
                    weight = w(i)*delta^2*pi;
                    maxw = max(maxw,weight);
                    sum = sum + weight*f1(r(i),a,delta)* ...
                        f2(theta(i),alpha(k,l),beta(k,l))*M;
                end
            elseif strcmp(discretization,'gauss-legendre')
                for i=1:n1
                    for j=1:n2
                        r = z(i)*delta;
                        theta = 2*pi*eta(j);
                        weight = w1(i)*w2(j)*2*pi*delta*r;
                        maxw = max(maxw,weight);
                        sum = sum + weight*f1(r,a,delta)* ...
                            f2(theta,alpha(k,l),beta(k,l))*M;
                    end
                end
            end
            Tbar_all(k,l) = sum;
        end
    end
end
Tbar = max(Tbar_all,[],'all');
tSpatial_initial = toc(tStart); % measure time for initial spatial discretization
    
if strcmp(tbound,'pessimistic')
    % Theoretical (pessimistic) bound
    tau = sspcoef * min(1/b, 1/(maxw*k1*k2*n*M+c));
elseif strcmp(tbound,'improved') || strcmp(tbound,'adaptive')
    % new time-step estimate
    tau = sspcoef * min(1/b, 1/(Tbar+c));
elseif strcmp(tbound,'given')
    if nargin > 8
        tau = tau_fixed;
    end
end
timestep = tau;

%% Time stepping procedure
t = 0;
numstep = 1;
finalstep = 0;
W = 0;
temp = 10^3; % choosing a big initial value to calculate the minimum
             % adaptive time step

% plot initial conditions
if strcmp(output{1},'plot_fig') || strcmp(output{1},'plot_error')
    fig = figure();
    fig.Renderer = 'painters';
    plotting(X,Y,S(1+gh:N-gh,1+gh:N-gh),I(1+gh:N-gh,1+gh:N-gh), ...
        R(1+gh:N-gh,1+gh:N-gh),t,S0,I0,S0+I0,tol,output)
    pause(0.01)
else
    fig = NaN;
end

tSpatial = 0;
tMethod = 0;
tStart_method = tic; % start clock
while t < tmax
    
    % always finish at the tmax time
    if t + tau >= tmax
        tau = tmax - t;
        finalstep = 1;
    end
    
    % The numerical method
    if strcmp(discretization,'num_integral')
        %% spatial discretization
        tStart_spatial = tic; % start clock
        
        Tinterior = computeT(XX,YY,I,X,Y,delta,a,alpha,beta,gridpts,N, ...
            gh,maxw,quad_pts,interpol,'gauss-legendre');
        
        tSpatial = tSpatial + toc(tStart_spatial); % measure time
                                           
        %% temporal discretization
        Snew = S(1+gh:N-gh,1+gh:N-gh).*exp(-tau*Tinterior-c*tau);
        %Snew=S(1+gh:N-gh,1+gh:N-gh).*(1 - tau*(Tinterior+c)+tau^2*(Tinterior+c).*(Tinterior+c)/2);
        Rnew = R(1+gh:N-gh,1+gh:N-gh)+b*tau*I(1+gh:N-gh,1+gh:N-gh) + ...
            tau*c*Snew;
        %Rnew=R(1+gh:N-gh,1+gh:N-gh)+b*tau*I(1+gh:N-gh,1+gh:N-gh)+tau*c*S(1+gh:N-gh,1+gh:N-gh);
        %Inew=(S(1+gh:N-gh,1+gh:N-gh)+I(1+gh:N-gh,1+gh:N-gh)+R(1+gh:N-gh,1+gh:N-gh))-Snew-Rnew;
        Inew = initsum - Snew - Rnew;
    else        
        % preallocation of new solution
        Z = zeros(3,gridpts,gridpts,s+1);
        
        % previous solution
        Y0(1,:,:) = S(1+gh:N-gh,1+gh:N-gh);
        Y0(2,:,:) = I(1+gh:N-gh,1+gh:N-gh);
        Y0(3,:,:) = R(1+gh:N-gh,1+gh:N-gh);
                                             
        % compute adaptive step size (only when FE is used)
        if strcmp(tbound,'adaptive') && strcmp(method,'FE')
            tStart_spatial = tic; % start clock for spatial discretization
            currentT(1+gh:N-gh,1+gh:N-gh) = computeT(XX,YY,I,X,Y, ...
                delta,a,alpha,beta,gridpts,N,gh,maxw,quad_pts,interpol, ...
                discretization);
            tSpatial = tSpatial + toc(tStart_spatial); % stop clock
            if finalstep == 0
                tau = min(1/b, 1/(max(max(currentT))+c));
                timestep = min(temp,tau);
                temp = timestep;                
            end
        end
        
        % compute stage values
        for ii=1:s+1
            Z(:,:,:,ii) = RK_delta(ii)*Y0(:,:,:);
            for jj = 1:ii-1
                if RK_beta(ii,jj) > 1e-14
                    II = zeros(N);
                    II(1+gh:N-gh,1+gh:N-gh) = Z(2,:,:,jj);
                    tStart_spatial = tic; % start clock for spatial discretization
                    Tinterior=computeT(XX,YY,II,X,Y,delta,a,alpha,beta, ...
                       gridpts,N,gh,maxw,quad_pts,interpol,discretization);
                    tSpatial = tSpatial + toc(tStart_spatial); % stop clock
                    Z(:,:,:,ii) = Z(:,:,:,ii) + ...
                        RK_beta(ii,jj)*(sspcoef*Z(:,:,:,jj) + ...
                        tau*rhs_fun(Z(:,:,:,jj),Tinterior,b,c,gridpts, ...
                        gridpts));
                end
            end
        end
        
        % update solution
        Snew(:,:) = Z(1,:,:,end);
        Inew(:,:) = Z(2,:,:,end);
        Rnew(:,:) = Z(3,:,:,end);
    end
    
    tMethod = tMethod + toc(tStart_method); % stop clock for time-stepping method
    
    % sanity checks
    if strcmp(property_checks,'checks_on')
        Wnew = errorcheck(gridpts,gridpts,Snew,Inew,Rnew, ...
            S(1+gh:N-gh,1+gh:N-gh),I(1+gh:N-gh,1+gh:N-gh), ...
            R(1+gh:N-gh,1+gh:N-gh),numstep,tol);
        W = max(W,Wnew);
    end
    
    % update
    tStart_method = tic; % restart clock
    S(1+gh:N-gh,1+gh:N-gh) = Snew(:,:);
    I(1+gh:N-gh,1+gh:N-gh) = Inew(:,:);
    R(1+gh:N-gh,1+gh:N-gh) = Rnew(:,:);
    t = t + tau;
    tMethod = tMethod + toc(tStart_method); % stop clock for time-stepping method
    
    % Plotting
    if strcmp(output{1},'plot_fig') || strcmp(output{1},'plot_error')
        plotting(X,Y,S(1+gh:N-gh,1+gh:N-gh),I(1+gh:N-gh,1+gh:N-gh), ...
            R(1+gh:N-gh,1+gh:N-gh),t,S0,I0,S0+I0,tol,output)
        pause(0.01)
    end
    
    tStart_method = tic; % restart clock
    numstep = numstep + 1;

end
tMethod = tMethod + toc(tStart_method); % stop clock for time-stepping method

if strcmp(property_checks,'checks_on')
    if W==0
        disp('The solution satisfies the qualitative conditions!')
    else
        disp('The solution violates the qualitative conditions!')
    end
end

% return values w.r.t physical domain
S = S(1+gh:N-gh,1+gh:N-gh);
I = I(1+gh:N-gh,1+gh:N-gh);
R = R(1+gh:N-gh,1+gh:N-gh);

CPUtimes = [tSpatial_initial, tSpatial, tMethod];

end

function plotting(X,Y,S,I,R,t,Smax,Imax,Rmax,tol,plot)
for i=1:length(S(:,1))
    for j=1:length(S(1,:))
        if S(i,j)<-tol
            S(i,j)=-1;
        end
        if I(i,j)<-tol
            I(i,j)=-1;
        end
        if R(i,j)<-tol
            R(i,j)=-1;
        end
    end
end

if strcmp(plot,'plot_fig')
    subplot(1,3,1)
    surf(X,Y,S)
    set(gca,'FontSize',12)
    if nargin>7
        axis([0,1,0,1,0,Smax])
    end
    subplot(1,3,2)
    surf(X,Y,I)
    set(gca,'FontSize',12)
    title(['$$t = \,$$',num2str(t)],'FontSize',16,'Interpreter','latex')
    if nargin>7
        axis([0,1,0,1,0,Imax])
    end
    subplot(1,3,3)
    surf(X,Y,R)
    set(gca,'FontSize',12)
    if nargin>7
        axis([0,1,0,1,0,Rmax])
    end
    set(gca,'DefaultAxesFontSize',12)
    
    colormap('default')
end

if strcmp(plot,'plot_error')
    surf(X,Y,S)
    hc = get(gcf,'children'); set(hc, 'fontsize', 12);
    title(['$$t = \,$$',num2str(t)],'FontSize',16,'Interpreter','latex')
    %if nargin>7
    axis([0,1,0,1,0,Smax+Imax+Rmax])
    %end
    axis('square')
    view(2)
end

end

function W = errorcheck(H,K,S,I,R,Sr,Ir,Rr,j,tol)
W = 0;
for i=1:H
    for l=1:K
        if S(i,l)<-tol
            disp('Negative value in S!')
            i
            l
            j
            S(i,l)
            W=1;
        end
        if I(i,l)<-tol
            disp('Negative value in I!')
            i
            l
            j
            I(i,l)
            W=1;
        end
        if R(i,l)<-tol
            disp('Negative value in R!')
            i
            l
            j
            R(i,l)
            W=1;
        end
        if S(i,l)>Sr(i,l)+tol
            disp('S is increasing!')
            i
            l
            j
            S(i,l)-Sr(i,l)
            W=1;
        end
        if R(i,l)<Rr(i,l)-tol
            disp('R is decreasing!')
            R(i,l)-Rr(i,l)
            i
            l
            j
            W=1;
        end
        if abs(S(i,l)+I(i,l)+R(i,l) - (Sr(i,l)+Ir(i,l)+Rr(i,l))) > tol
            disp('Conservation violated')
            S(i,l)+I(i,l)+R(i,l) - (Sr(i,l)+Ir(i,l)+Rr(i,l))
            i
            l
            j
            W=1;
        end
    end
end
% if W == -1
%     W = 0;
% end

end

function [out]=outhandle(x,y,A,B)
% Handles the cases in which the point (x,y) lies outside of our domain
% "out" is the output, the coordinate which should be considered, or gives back 0 in the no boundary case.
% pr - if '1', then it runs the programm for the Neumann boundary case
%      if '2', then it runs the no boundary case
% if pr==1
%     if x<=A && y<=B && 0<=x && 0<=y
%         out=[x,y];
%     end
%     if x<0
%         if y<0
%             out=[0,0];
%         end
%         if 0<=y && y<=B
%            out=[0,y];
%         end
%         if y>B
%            out=[0,B];
%         end
%     end
%     if 0<=x && x<=A
%         if y<0
%             out=[x,0];
%         end
%         if y>B
%            out=[x,B];
%         end
%     end
%     if x>A
%         if y<0
%             out=[A,0];
%         end
%         if 0<=y && y<=B
%            out=[A,y];
%         end
%         if y>B
%            out=[A,B];
%         end
%     end
% end

% if pr==2
if x<=A && y<=B && 0<=x && 0<=y
    out=[x,y];
end
if x<0
    if y<0
        out=0;
        %disp('--')
    end
    if 0<=y && y<=B
        out=0;
        %disp('-0')
    end
    if y>B
        out=0;
        %disp('-+')
    end
end
if 0<=x && x<=A
    if y<0
        out=0;
        %disp('0-')
        %             disp('Here.')
    end
    if y>B
        out=0;
        %disp('0+')
    end
end
if x>A
    if y<0
        out=0;
        %disp('+-')
    end
    if 0<=y && y<=B
        out=0;
        %disp('+0')
    end
    if y>B
        out=0;
        %disp('++')
    end
end
% end

end

function val = f1(r, a, delta)
% Determines the values of f_1.
val = a*(-r+delta);

% if r>delta
%     val=0;
% end
end

function val = f2(theta, alpha, beta)
% Determines the values of f_2.
val = beta*(sin(pi/2 + theta - alpha) + 11/10);
end

function F = rhs_fun(Y,T,b,c,H,K)
% Computes the RHS for the ODE system of the SIR problem
%   S'= -S*T - c*S;;
%   I' = S*T - b*I;
%   R' = b*I + c*S;

F = zeros(3,H,K);
S(:,:) = Y(1,:,:);
I(:,:) = Y(2,:,:);

F(1,:,:) = -S.*T - c*S;
F(2,:,:) =  S.*T - b*I;
F(3,:,:) =  b*I + c*S;
end

function Tinterior = computeT(XX,YY,I,X,Y,delta,a,alpha,beta,gridpts,N, ...
    gh,maxw,quad_pts,interpol,discretization)
%    gh,maxw,pr2,pr3,type,interpol_type)                  
T = zeros(N);

if strcmp(discretization,'4point')
    % The 4 point method
    
    Iinterp1 = interp2(XX,YY,I,X+1/2*delta,Y+1/2*delta,interpol);
    Iinterp2 = interp2(XX,YY,I,X-1/2*delta,Y+1/2*delta,interpol);
    Iinterp3 = interp2(XX,YY,I,X-1/2*delta,Y-1/2*delta,interpol);
    Iinterp4 = interp2(XX,YY,I,X+1/2*delta,Y-1/2*delta,interpol);

    for k = 1+gh:N-gh
        for l = 1+gh:N-gh
            T(k,l) = maxw*f1(delta/sqrt(2),a,delta)* ...
                (Iinterp1(k-gh,l-gh)*f2(pi/4,alpha(k,l),beta(k,l))+ ...
                Iinterp2(k-gh,l-gh)*f2(3*pi/4,alpha(k,l),beta(k,l))+ ...
                Iinterp3(k-gh,l-gh)*f2(5*pi/4,alpha(k,l),beta(k,l))+ ...
                Iinterp4(k-gh,l-gh)*f2(7*pi/4,alpha(k,l),beta(k,l)));
        end
    end
else
    val = zeros(gridpts,gridpts);
    
    n1 = quad_pts;
    n2 = 2*n1;
    n = n1*n2;
    % get quadrature points
    if strcmp(discretization,'elhay-kautsky')
        [w, ~, ~, r, theta] = disk_rule(n1, n2, 0, 0, delta);
    elseif strcmp(discretization,'gauss-legendre')
        [z, w1] = gauss_legendre_rule_compute(n1, 0, 1);
        [eta, w2] = gauss_legendre_rule_compute(n2, 0, 1);
    end
    
    xpoints = X;
    ypoints = Y;
    
    for k=1:gridpts
        for l=1:gridpts
            % The values which will be written inside the integral
            %The coordinates of the point (i,l):
            pointx = xpoints(k,l);
            pointy = ypoints(k,l);
            
            % Let j be the number of the circle, and k is the number of the
            % given point on this circle.
            if strcmp(discretization,'elhay-kautsky')
                %The coordinates of the n1 x n2 points
                xcoord=zeros(n,1);
                ycoord=zeros(n,1);
                
                for i=1:n
                    xplus = r(i)*cos(theta(i));
                    yplus = r(i)*sin(theta(i));
                    
                    xcoord(i) = pointx + xplus;
                    ycoord(i) = pointy + yplus;
                end
            elseif strcmp(discretization,'gauss-legendre')
                %The coordinates of the n1 x n2 points
                xcoord=zeros(n1,n2);
                ycoord=zeros(n1,n2);
                
                for i=1:n1
                    for j=1:n2
                        % The value which will be added to the coordinate of point
                        % (i,l)
                        
                        xplus = z(i)*delta*cos(2*pi*eta(j));
                        yplus = z(i)*delta*sin(2*pi*eta(j));
                        
                        xcoord(i,j) = pointx + xplus;
                        ycoord(i,j) = pointy + yplus;
                        
                    end
                end
                
            end
            
            interpolating_whole_function = 0;
            if interpolating_whole_function == 1
                % interpolating the whole function
                V = zeros(N);
                for ind1=1:gridpts
                    for ind2=1:gridpts
                        radius=sqrt((pointx - xpoints(ind1,ind2))^2+(pointy - ypoints(ind1,ind2))^2);
                        thetta=atan2((pointy - ypoints(ind1,ind2)),(pointx - xpoints(ind1,ind2)));
                        
                        V(ind1+gh,ind2+gh)=radius*(a*(-radius*delta + delta))*(beta(k,l)*(sin(thetta+alpha(k,l))+1))*I(ind1,ind2);
                    end
                end
                Ip=interp2(XX,YY,V,xcoord,ycoord,interpol);
            else
                Ip=interp2(XX,YY,I,xcoord,ycoord,interpol);
            end
            
            %The construction:
            sum = 0;
            if strcmp(discretization,'elhay-kautsky')
                Ip = reshape(Ip,[],n);
                for i=1:n
                    weight = w(i)*delta^2*pi;
                    sum = sum + weight*f1(r(i),a,delta)* ...
                        f2(theta(i),alpha(k,l),beta(k,l))*Ip(i);
                end
            elseif strcmp(discretization,'gauss-legendre')
                for i=1:n1
                    for j=1:n2
                        weight = w1(i)*w2(j)*2*pi*delta^2*z(i);
                        sum = sum + weight*f1(z(i)*delta,a,delta)* ...
                            f2(2*pi*eta(j),alpha(k,l),beta(k,l))* ...
                            Ip(i,j);
                    end
                end
            elseif interpolating_whole_function == 1
                Ip = reshape(Ip,[],n);
                for i=1:n
                    weight = w(i)*delta^2*pi;
                    sum = sum + weight*Ip(i);
                end
            end
            
            val(k,l) = sum;
        end
    end
    
    T(1+gh:N-gh,1+gh:N-gh) = val;
end

Tinterior = T(1+gh:N-gh,1+gh:N-gh);

end
