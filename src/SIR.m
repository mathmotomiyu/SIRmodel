% file SIR.m
% ===========
% PART 1
% ===========
clc
clear
close all

L = 60;
S0 = 0.999; I0 = 0.001; R0 = 0;
x0 = [S0;I0;R0]; % initial  condition
params.u1 = 0; % lower bound for control
params.u2 = 0.1; % upper bound for control: constant M
params.m = 0; % constant in control system
params.d = 0.1; % constant in control system
params.R0 = 10; % constant in control system: constant cN/d
params.c = params.R0 * params.d;
params.a = 0.2; % constant in cost

h = 0.0001;
params.h = h; % time grid step
t = 0:h:L; % time grid, row vector
n = length(t);
x = zeros(3,n); % state vector - vector of 3 variables S, I, R
x1 = zeros(3,n); % state vector to be used for ρ loop (Armijo method)

p = zeros(3,n); % adjoint state vector
unew = zeros(1,n); % control vector corresponding to u(k+1)
u = zeros(1,n); % control vector to be used for ρ loop (Armijo method)
u0 = params.u2/1.0; % initial control value
uold = u0*ones(1,n); % control vector corresponding to u(k)

eps = 0.00001; % precision ε for the gradient algorithm
maxit = 500; % max. no. iterations – gradient algorithm disp('RO data')
roin = 1; % initial value for gradient steplength ρ
bro = 0.6; % b parameter for ρ loop
eps1 = 0.00001; % precision ε1 for steplength ρ
maxro = 20; % max. no. of iterations – ρ loop
ro = roin; % initialization of steplength ρ
flag1 = 0; % convergence indicator for gradient algorithm

grad = zeros(1,maxit);
cost = zeros(1,maxit);

% ===========
% PART 2
% ===========
% gradient loop starts
for iter = 1:maxit
    fprintf('iter = %d \n',iter);
    if iter == 1
        % S1 : solve the state equation by Runge–Kutta method
        x(:,1) = x0; % the initial condition
        for i = 1:n-1
            x(:,i+1) = RK41(params,t(i),x(:,i),uold(i));
        end
        disp('State Problem solved');
        
        Q = x(2,:)+params.a*uold.*uold;
        cvold = trapz(t,Q);
        jj = 1;
        cost(jj) = cvold;
    end

    % ===========
    % PART 3
    % ===========
    % S2 : solve adjoint equation
    p(:,n) = 0;
    for j = 1:n-1
        i = n-j;
        p(:,i) = RK43(params,t(i+1),x(:,i+1),uold(i+1),p(:,i+1));
    end
    disp('AE solved');
    % S3 : compute the gradient
    w = x(1,:) .* (p(1,:) - p(3,:)) + 2 * params.a * uold;
    disp('GRADIENT COMPUTED');

    % ===========
    % PART 4
    % ===========
    normg = sqrt(sum(w.^2)/n); % compute the gradient norm
    grad(jj) = normg;
    % verify SC1
    if normg < eps
        disp('CONVERGENCE by GRADIENT')
        flag1 = 1;
        break
    end

    % ===========
    % PART 5
    % ===========
    % S4 : FIT RO
    robar = ro;
    flag2 = 0;
    flag3 = 0;
    % start loop to fit ro
    for count = 1:maxro
        for i = 1:n
            temp = uold(i) - robar*w(i);
            u(i) = Proj(params,temp);
        end
        % solve state equation for input u and get state x1
        x1(:,1) = x0;
        for j = 1:n-1
            x1(:,j+1) = RK41(params,t(j),x1(:,j),u(j));
            %x1(j+1) = x1(j) + h*F1(params,t(j),x1(j),u(j));
        end
        % test for SC4
        if robar < eps1
            flag3 = 1;
            break % leave loop for count . . .
        end
        % compute current cost value and test ¥bar{ρ}
        Q1 = x1(2,:)+params.a*u.*u;
        cv = trapz(t,Q1);
        if cv >= cvold % no decrease of cost for minimization
            robar = bro * robar;
        else
            flag2 = 1;
            break % leave loop for count . . .
        end
    end % for count
    fprintf('counts = %d \n',count);
    
    % ===========
    % PART 6
    % ===========
    if flag3 == 1
        disp('CONVERGENCE BY RO')
        break % leave loop for iter . . . (SC4 satisfied)
    end
    if flag2 == 1
        cvnew = cv;
        jj = jj + 1;
        cost(jj) = cvnew;
        unew = u;
        ro = robar;
        % verify SC2
        if abs(cvnew - cvold) < eps
            disp('CONVERGENCE by COST')
            flag1 = 1;
            break % leave the loop for iter . . . (SC2 satisfied)
        end
        % verify SC3
        d = uold - unew;
        dif = sqrt(sum(d.^2));
        if dif < eps
            disp('CONVERGENCE by CONTROL')
            flag1 = 1;
            break % leave the loop for iter . . . (SC3 satisfied)
        end
    else         
        error('NO CONVERGENCE FOR RO') % STOP PROGRAM
    end

    % ===========
    % PART 7
    % ===========
    % prepare a new iteration
    uold = unew;
    cvold = cvnew;
    x = x1;
end % for iter
if (flag1 == 1) || (flag3 == 1)
    figure(1)
    %subplot(2,2,1)
    plot(t,u,'LineWidth',2); grid
    hold on
    xlabel('\bf t','FontSize',16)
    p(:,n) = 0;
    for j = 1:n-1
        i = n-j;
        p(:,i) = RK43(params,t(i+1),x(:,i+1),uold(i),p(:,i+1));

    end
    graph = (x(1,:) .* (p(3,:) - p(1,:))) / (2*params.a);
    plot(t,graph,'LineWidth',2); grid
    ylim([-0.01,2*params.u2]);
    legend('u(t)','S(t)(p3(t)-p1(t))/2a');

    figure(2)
    %subplot(2,2,2)
    plot(t,x1(1,:),'LineWidth',2); grid
    hold on
    plot(t,x1(2,:),'LineWidth',2);
    plot(t,x1(3,:),'LineWidth',2);
    xlabel('\bf t','FontSize',16)
    ylabel('\bf x(t)','FontSize',16)
    legend('S(t)','I(t)','R(t)');

    figure(3)
    %subplot(2,2,3)
    plot(cost(1:iter),'LineWidth',2); grid
    legend('cost');

    figure(4)
    %subplot(2,2,4)
    plot(grad(1:iter),'LineWidth',2); grid
    legend('gradient size');
else
    error('NO CONVERGENCE FOR DESCENT METHOD')
end
save grad.txt grad -ascii % save vector grad into file grad.txt
save cost.txt cost -ascii % save vector cost into file cost.txt
disp('END OF JOB')
