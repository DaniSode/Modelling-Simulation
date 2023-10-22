clear all; close all; clc
%% To note is that we have used parts of the code structure from PSS 6. 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                               %
%                          Assignment 3                         %                
%                   Modelling and simulation                    %
%        Written by Johannes Lundahl and Daniel Söderqvist      %
%                        __ oktober 2023                        %
%                                                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%% Task 2

%% Defining functions for RK1, RK2, RK4 which are tested and ploted 
% Defining initial values
lambda = -2;
x0 = 1;

% Simulation parameters:
tFinal  = 2;
deltaT   = 0.1;

% The exact Solution:
tExact = 0:10^-6:tFinal;
xExact = x0*exp(lambda*tExact);


% RK1

% Butchers table
c1= 0;
b1 = 1;
BT_RK1 = [c1 0;
          0 b1];

% Testing RK1
[tRK1,xRK1] = RK1(tFinal, deltaT, x0, lambda, BT_RK1);


% RK2

%Butchers table
a = 1/2;
c1= 0;
c2 = 1/2;
b1 = 0;
b2 = 1;    
BT_RK2 = [c1 0 0;
      c2 a 0;
      0 b1 b2];

% Testing RK2
[tRK2,xRK2] = RK2(tFinal, deltaT, x0, lambda, BT_RK2);


% RK4

% Butchers table
a1 = 1/2;
a2 = 1/2;
a3 = 1;
c1= 0;
c2 = 1/2;
c3 = 1/2;
c4 = 1;
b1 = 1/6;
b2 = 1/3;
b3 = 1/3;
b4 = 1/6;
BT_RK4 = [c1 0 0 0 0;
          c2 a1 0 0 0;
          c3 0 a2 0 0;
          c4 0 0 a3 0;
          0 b1 b2 b3 b4];

% Testing RK4
[tRK4,xRK4] = RK4(tFinal, deltaT, x0, lambda, BT_RK4);

% Plotting:
figure
subplot(1,3,1)
plot(tExact, xExact,'markersize',20)
hold on
plot(tRK1,xRK1,'marker','.','markersize',10)
xlabel('t');
ylabel('x');
title('RK1')
legend('Exact','RK1')

subplot(1,3,2)
plot(tExact, xExact,'markersize',20)
hold on
plot(tRK2,xRK2,'marker','.','markersize',10)
xlabel('t');
ylabel('x');
title('RK2')
legend('Exact','RK2')

subplot(1,3,3)
plot(tExact, xExact,'markersize',20)
hold on
plot(tRK4,xRK4,'marker','.','markersize',10)
xlabel('t');
ylabel('x');
title('RK4')
legend('Exact','RK4')


%% Compute error with different deltaT and plot
deltaTs = logspace(-3,-1,3);
errorRK1 = zeros(1, length(deltaTs));
errorRK2 = zeros(1, length(deltaTs));
errorRK4 = zeros(1, length(deltaTs));

for i = 1:length(deltaTs)
    [t1,xRK1] = RK1(tFinal, deltaTs(i), x0, lambda, BT_RK1);
    [t2,xRK2] = RK2(tFinal, deltaTs(i), x0, lambda, BT_RK2);
    [t4,xRK4] = RK4(tFinal, deltaTs(i), x0, lambda, BT_RK4);
    xExact_1 = x0*exp(lambda*t1);
    xExact_2 = x0*exp(lambda*t2);
    xExact_4 = x0*exp(lambda*t4);
    errorRK1(i) = vpa(norm(xRK1(end)-xExact_1(end),inf),20);
    errorRK2(i) = vpa(norm(xRK2(end)-xExact_2(end),inf),20);
    errorRK4(i) = vpa(norm(xRK4(end)-xExact_4(end),inf),20);
end

figure
loglog(deltaTs,errorRK1,'marker','.','markersize',15)
hold on; grid on
loglog(deltaTs,errorRK2,'marker','.','markersize',15)
loglog(deltaTs,errorRK4,'marker','.','markersize',15)

set(gca,'XDir','reverse')
legend('RK1', 'RK2', 'RK4')
xlabel('dt')
ylabel('Error')
title('Error vs different delta Ts')


%% Stability with print and plot

% Define syms and needed variables
orders = [1 2 4];
solution = zeros(length(orders));
syms lambda

% Calculate sum and look for stability using formula from 7.41 page 170
for i = 1:length(orders)
    summation = 0; 
    for j = 0:orders(i)   
        summation = summation + ((lambda*deltaT)^j/factorial(j));
    end
    summation = abs(summation);
solution(i) = double(solve(summation==1, lambda<0));
end
solution = solution(:,1);
fprintf('The minimum lambda for which the solution is still stable\nis shown down below, for RK1, RK2 and RK4 respective:\n'); disp(solution)

% Showing solution as graphs
maxl = -50;
lam = linspace(maxl,0,100);
summation = zeros(length(lam),length(orders));

for i = 1:length(orders)
    for l = 1:length(lam) 
        for k=0:orders(i)
            summation(l,i) = summation(l,i) + ((lam(l)*deltaT(end))^k)/factorial(k);
        end
        summation(l,i) = abs(summation(l,i));
    end
end

figure
plot(lam,summation(:,1))
hold on; grid on
plot(lam,summation(:,2))
plot(lam,summation(:,3))
plot(lam, ones(length(lam)),'k')
scatter(solution, ones(length(solution)), 'ko')
xlabel('lambda'); ylabel('S'); title('Stability graphs') 
legend('RK1', 'RK2', 'RK4', 'Threshold')


%% Solving van der pol with ode45
y0 = [0 1];
tf = [0 25];
[t,y] = ode45(@vanderpol, tf, y0, odeset());

figure
subplot(2,1,1)
plot(t,y)
grid on
title('Solution of van der Pol Equation with ODE45');
xlabel('t');
ylabel('y');
legend('y_1','y_2')

subplot(2,1,2)
plot(diff(t))
grid on; axis tight
title('Discrete time-step of ode45 function');
xlabel('t');
ylabel('deltaT');
legend('Discrete time-steps')


%% Solving van der pol with RK4

% Butchers table
a1 = 1/2;
a2 = 1/2;
a3 = 1;
c1= 0;
c2 = 1/2;
c3 = 1/2;
c4 = 1;
b1 = 1/6;
b2 = 1/3;
b3 = 1/3;
b4 = 1/6;
BT_RK4 = [c1 0 0 0 0;
          c2 a1 0 0 0;
          c3 0 a2 0 0;
          c4 0 0 a3 0;
          0 b1 b2 b3 b4];

% Solving
% Try different delta Ts
deltaTt = 0.01;
[tRK4,xRK4] = RK4vanderpol(tf(2), deltaTt, y0, BT_RK4);
   
% Plot
figure
subplot(2,1,1)
[t,y] = ode45(@vanderpol, tf, y0, odeset('AbsTol',1e-8,'RelTol',1e-8));
plot(tRK4,xRK4,'b',t,y,'r')
grid on; axis tight
title('Solution of van der Pol Equation vs RK4');
xlabel('t');
ylabel('y');
legend('xRK4_1','xRK4_2', 'y_1', 'y_2')

subplot(2,1,2)
plot(diff(t))
grid on; axis tight
title('Discrete timestep of ode45 with higher tolerance');
xlabel('t');
ylabel('dt');
legend('Discrete timestep')



%% IRK

% Butchers table
BT_IRK4 = [1/2 - sqrt(3)/6, 1/4, 1/4 - sqrt(3)/6;
           1/2 + sqrt(3)/6, 1/4 + sqrt(3)/6, 1/4;
           0,               1/2,             1/2];

% Define needed syms
syms lambda t deltaT xk

% Define r(K,xk,u)
s = 2; %order dim
n = 1;  %state space dim
K = sym('K',[s,n]);
A = sym('A', [s,s]);

% Define r and dr 
r = [testfunction(xk + deltaT*A(1,1)*K(1,:) + deltaT*A(1,2)*K(2,:),lambda)-K(1,:);
     testfunction(xk + deltaT*A(2,1)*K(1,:) + deltaT*A(2,2)*K(2,:), lambda)-K(2,:)];

dr = jacobian(r,K);

% Define function
matlabFunction(r,dr, 'file', 'rdrIRK','vars',{deltaT,xk,K,A,lambda});

% Define needed variables
lambda = -2;
deltaT= 0.1;
tf = 10;
xk = 1;
A = BT_IRK4(1:end-1,2:end);
b = BT_IRK4(3,2:end);
nIRK = tf / deltaT;
tIRK = zeros(nIRK,n);
xIRK = zeros(nIRK,n);
xIRK(1,:) = xk;
tol = 10^-6;
alfa = 1;

% Using formula on page 180 to calcualte r until tolerance reached 
for j = 1:nIRK
    tIRK(j+1,1) = tIRK(j) + deltaT;
    K_j  = [xIRK(j,1);xIRK(j,1)];
    [r,dr] = rdrIRK(deltaT,xIRK(j,1),K_j,A,lambda);
  while abs(r) > tol
      [r,dr] = rdrIRK(deltaT,xIRK(j,1),K_j,A,lambda);
      deltaK = -dr\r;
      K_j  = K_j + alfa*deltaK; 
  end
  xIRK(j+1,:) = xIRK(j,:) + deltaT*K_j(1)*b(1) + deltaT*K_j(2)*b(2); 
end

% Define needed syms
syms lambda t deltaT real
xk = sym('xk',[2,1],'real');
s = 2; %order dim
n = 2;  %state space dim
K = sym('K',[s,n],'real');
A = sym('A', [s,s],'real');

% Define r and dr for vdp
r_vdp = [vanderpol(t,xk + deltaT*A(1,1)*K(:,1) + deltaT*A(1,2)*K(:,2))-K(:,1);
        vanderpol(t, xk + deltaT*A(2,1)*K(:,1) + deltaT*A(2,2)*K(:,2))-K(:,2)];
K = reshape(K,[numel(K),1]);
dr_vdp = jacobian(r_vdp,K);

% Define function
matlabFunction(r_vdp,dr_vdp, 'file', 'rdrIRKvdp','vars',{deltaT,xk,K,A,t});

% Define needed variables
deltaT2 = 0.01;
tf = 25;
xk = [1 0];
A = BT_IRK4(1:end-1,2:end);
b = BT_IRK4(3,2:end);
nIRK_vdp = tf / deltaT2;
tIRK_vdp = zeros(nIRK_vdp,n);
xIRK_vdp = zeros(n,nIRK_vdp);
xIRK_vdp(:,1) = xk;
tol = 10^-6;
alfa = 1;

% Using fomrula on page 180
for j = 1:nIRK_vdp
    tIRK_vdp(j+1,1) = tIRK_vdp(j) + deltaT2;
    K_j  = [xIRK_vdp(:,j)',xIRK_vdp(:,j)']';
    [r,dr] = rdrIRKvdp(deltaT2,xIRK_vdp(:,j),K_j, A, tIRK_vdp(j,1));
  while abs(r) > tol
      [r,dr] = rdrIRKvdp(deltaT2,xIRK_vdp(:,j),K_j, A, tIRK_vdp(j,1));
      deltaK = -dr\r;
      K_j  = K_j + alfa*deltaK;  
  end
    K_j = reshape(K_j,[],2);
    xIRK_vdp(:,j+1) = xIRK_vdp(:,j) + deltaT2*K_j(:,1)*b(1) + deltaT2*K_j(:,2)*b(2); 
end


%% Stability
maxl = -10000;
lam = linspace(maxl,0,100);
R = zeros(length(lam),1);
for l = 1:length(lam) 
        R_k = 1 + lam(l)*deltaT2*b*pinv(eye(2)-lam(l)*deltaT2*A)*ones(2,1);
        R(l) = norm(R_k);
end

figure
plot(lam, R)
xlabel('lambda'); ylabel('R'); title('Stability for IRK')

%% VS RK4

[tRK4_4,xRK4_4] = RK4vanderpol(tf, deltaT2, xk, BT_RK4);
[t,y] = ode45(@vanderpol, [0 tf], xk, odeset('AbsTol',1e-8,'RelTol',1e-8));

figure
subplot(2,1,1)
plot(tRK4_4,xRK4_4,'r')
hold on; grid on
plot(tIRK_vdp(:,1),xIRK_vdp,'b')
plot(t, y, 'k')
xlabel('Time t');
ylabel('Solution y');
legend('y_1 IRK4','y_2 IRK4','y_1 RK4','y_2 RK4', 'y_1 true', 'y_2 true')
title('Solution of van der Pol Equation with IRK4 and RK4');
axis tight

subplot(2,1,2)
plot(tRK4_4, abs(xRK4_4'-xIRK_vdp))
hold on; grid on
xlabel('t');
ylabel('y');
legend('Error y_1', 'Error y_2')
title('Error between IRK4 and ERK4');
axis tight






