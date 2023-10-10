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

%% Defining initial values
%unstable RK1 & RK2 = ||-20||
%unstable RK4 = ||-28||
lambda = -10;
x0 = 1;
%% Simulation parameters:
tFinal  = 5;
deltaT   = [0.001, 0.01, 0.1];

%% The exact Solution:
tExact = 0:0.01:tFinal;
xExact = x0*exp(lambda*tExact);

%% RK1
c1= 0;
b1 = 1;
%Butchers table
BT_RK1 = [c1 0;
          0 b1];
tRK1_alld = zeros(1,3);
xRK1_alld = zeros(1,3);
for i = 1:length(deltaT)
    [tRK1,xRK1] = RK1(tFinal, deltaT(i), x0, lambda, BT_RK1);
    tRK1_alld(1:length(tRK1),i) = tRK1;
    xRK1_alld(1:length(xRK1),i) = xRK1;
end

%% RK2
a = 1/2;
c1= 0;
c2 = 1/2;
b1 = 0;
b2 = 1;
      
BT_RK2 = [c1 0 0;
      c2 a 0;
      0 b1 b2];

tRK2_alld = zeros(1,3);
xRK2_alld = zeros(1,3);

for i = 1:length(deltaT)
    [tRK2,xRK2] = RK2(tFinal, deltaT(i), x0, lambda, BT_RK2);
    tRK2_alld(1:length(tRK2),i) = tRK2;
    xRK2_alld(1:length(xRK2),i) = xRK2;
end
%% RK4
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

%Butchers table
BT_RK4 = [c1 0 0 0 0;
          c2 a1 0 0 0;
          c3 0 a2 0 0;
          c4 0 0 a3 0;
          0 b1 b2 b3 b4];

tRK4_alld = zeros(1,3);
xRK4_alld = zeros(1,3);

for i = 1:length(deltaT)
    [tRK4,xRK4] = RK4(tFinal, deltaT(i), x0, lambda, BT_RK4);
    tRK4_alld(1:length(tRK4),i) = tRK4;
    xRK4_alld(1:length(xRK4),i) = xRK4;
end
%% Plotting:
figure
subplot(1,3,1)
plot(tExact, xExact,'markersize',20)
hold on
plot(tRK1,xRK1,'marker','.','markersize',10)
xlabel('t');
ylabel('x');
legend('Exact','RK1')

subplot(1,3,2)
plot(tExact, xExact,'markersize',20)
hold on
plot(tRK2,xRK2,'marker','.','markersize',10)
xlabel('t');
ylabel('x');
legend('Exact','RK2')

subplot(1,3,3)
plot(tExact, xExact,'markersize',20)
hold on
plot(tRK4,xRK4,'marker','.','markersize',10)
xlabel('t');
ylabel('x');
legend('Exact','RK4')

%% Compute error
err = zeros(3,length(deltaT));

err(1,3) = norm(xRK1_alld(21,3)-xExact(end),inf);
err(1,2) = norm(xRK1_alld(201,2)-xExact(end),inf);
err(1,1) = norm(xRK1_alld(end,1)-xExact(end),inf);
err(2,3) = norm(xRK2_alld(21,3)-xExact(end),inf);
err(2,2) = norm(xRK2_alld(201,2)-xExact(end),inf);
err(2,1) = norm(xRK2_alld(end,1)-xExact(end),inf);
err(3,3) = norm(xRK4_alld(21,3)-xExact(end),inf);
err(3,2) = norm(xRK4_alld(201,2)-xExact(end),inf);
err(3,1) = norm(xRK4_alld(end,1)-xExact(end),inf);

figure(2);
loglog(deltaT,err(1,:),'marker','.','markersize',15)
hold on
loglog(deltaT,err(2,:),'marker','.','markersize',15)
loglog(deltaT,err(3,:),'marker','.','markersize',15)

grid on
set(gca,'XDir','reverse')
legend('RK1', 'RK2', 'RK4')
xlabel('dt')
ylabel('Error')

%% Stability using R-formula
% syms lambda lambda_RK2 negative real
% my = lambda*deltaT(end);
% 
% b1 = BT_RK1(2,2:end);
% A = BT_RK1(1,2:end);
% R = 1 + my*b1.'*inv(eye(1)-my*A)*ones(1,1);
% lam1 = solve(R<=-1,lambda)
% 
% my_RK2 = lambda_RK2*deltaT(end);
% b2 = BT_RK2(3,2:end);
% A2 = BT_RK2(1:2,2:end);
% R2 = 1 + my_RK2*b2*((eye(2)-my_RK2*A2)^-1*ones(2,1));
% lam2 = solve(R<=-1,lambda_RK2)

%% Using S-function
orders = [1, 2, 4];
maxl = -50;
lam = linspace(maxl,0,100);
summation = zeros(length(lam),length(orders));

for i = 1:length(orders)
    for l = 1:length(lam) 
        for k=0:orders(i)
            summation(l,i) = summation(l,i) + ((lam(l)*deltaT(end))^k)/factorial(k);
        end
        summation(l,i) = abs(summation(l,i));
%         
%         if summation(l,i) < -1
%             fprintf('\nFor order %.f the solution becomes unstable at lambda = %.f\n\n' , orders(i), lam(l))
%             break
%         end
%         if lam == 0
%             disp('Increase max lambda')
%         end
    end
end

figure
plot(lam,summation(:,1))
hold on
plot(lam,summation(:,2))
plot(lam,summation(:,3))

%% Q3 - Solving van der pol with ode45
y0= [0 1];
tf = [0 25];
[t,y] = ode45(@vanderpol, tf,y0);

figure
subplot(2,1,1)
plot(t,y(:,1),'-o',t,y(:,2),'-o')
title('Solution of van der Pol Equation with ODE45');
xlabel('Time t');
ylabel('Solution y');
legend('y_1','y_2')

%% Solving vdp with RK4
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

%Butchers table
BT_RK4 = [c1 0 0 0 0;
          c2 a1 0 0 0;
          c3 0 a2 0 0;
          c4 0 0 a3 0;
          0 b1 b2 b3 b4];


[tRK4,xRK4] = RK4vanderpol(tf(2), deltaT(end), y0, BT_RK4);
   
hold on
subplot(2,1,2)
plot(tRK4,xRK4(:,1),'-o',t,y(:,2),'-o')
title('Solution of van der Pol Equation with RK4');
xlabel('Time t');
ylabel('Solution y');
legend('y_1','y_2')

%% Q4 -IRK

BT_IRK4 = [1/2 - sqrt(3)/6, 1/4, 1/4 - sqrt(3)/6;
           1/2 + sqrt(3)/6, 1/4 + sqrt(3)/6, 1/4;
           0,               1/2,             1/2];

syms lambda t deltaT xk

% Define r(K,xk,u)
s = 2; %order dim
n = 1  %state space dim
K = sym('K',[s,n]);
A = sym('A', [s,s]);

r = [f(xk + deltaT*A(1,1)*K(1) + deltaT*A(1,2)*K(2),u)-K(1);
     f(xk + deltaT*A(2,1)*K(1) + deltaT*A(2,2)*K(2), lambda)-K(2)];

dr = jacobian(r,K)

matlabFunction(r,dr, 'file', 'rdrIRK','vars',{deltaT,t,xk,K,A,lambda});
clear delaT xk r dr 


xdot = lamda *x

f(x,u,t)












