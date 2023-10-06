clear all; close all; clc
%% To note is that we have used parts of the code structure from PSS 6. 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                               %
%                          Assignment 2                         %                
%                   Modelling and simulation                    %
%        Written by Johannes Lundahl and Daniel Söderqvist      %
%                       27 september 2023                       %
%                                                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%% 1.3a)
% Defining variables and least squares function
a0 = -2.923;
b0 = 7.18;
N = 100;
my = 0;
sigma = 3.8;
u = unifrnd(0,50,[N 1]);
u = sort(u);
e = normrnd(my, sigma, [N 1]);
y = a0 + b0*u + e;



%% 1.3b)
% Defining variables and 2nd order polynomial function
c0 = 2.8;
my2 = 0;
sigma2 = 12.8;
eps = normrnd(my2, sigma2, [N 1]);
y_real = a0 + b0*u + c0*(u.^2) + eps;

% Defining theta and y_hat first case
phi = [ones(N, 1), u];
theta_hat = pinv(phi'*phi)*phi'*y;
y_hat = theta_hat'*phi';
fprintf('#####################\nParameters all cases\n#####################\n\n')
fprintf(['The parameters from first case:\n\n'...
    'a_0 = %.3f\nb_0 = %.3f\n\n'], theta_hat(1),theta_hat(2))

% Plotting real linear data vs least square fit
rng(1)
plot(u,y,'ok', 'MarkerSize', 5, 'MarkerFaceColor','b')
hold on 
plot(u, y_hat,'LineWidth',1.2 ,'color', 'r')

% Defining theta and y_hat second case
phi2 = [ones(N, 1), u, u.^2];
theta_hat2 = pinv(phi2'*phi2)*phi2'*y_real;
y_hat2 = theta_hat2'*phi2';
fprintf(['The parameters from secong case:\n\n'...
    'a_0 = %.3f\nb_0 = %.3f\nc_0 = %.3f\n\n'], theta_hat2(1), theta_hat2(2), theta_hat2(3))

% Plotting real poly data vs 2nd order poly fit function
figure
rng(1)
plot(u,y_real,'ok',  'MarkerSize', 5,'MarkerFaceColor','b')
hold on
plot(u, y_hat2, 'LineWidth',1.2 ,'color', 'r')

% Defining theta and y_hat third case
phi3 = [ones(N, 1), u];
theta_hat3 = pinv(phi3'*phi3)*phi3'*y_real;
y_hat3 = theta_hat3'*phi3';
fprintf(['The parameters from third case:\n\n'...
    'a_0 = %.3f\nb_0 = %.3f\n\n\n'], theta_hat3(1), theta_hat3(2))

% Plotting real poly data vs least square fit from first plot
figure
plot(u,y_real,'ok', 'MarkerSize', 5, 'MarkerFaceColor','b')
hold on
plot(u, y_hat3, 'LineWidth',1.2 ,'color', 'r')

% Calculating the total residuals and printing to screen
eps_curve = (1/N)*sum((y_real - y_hat2').^2);
eps_linear = (1/N)*sum((y_real - y_hat3').^2); 
fprintf('#####################\n Resulting residuals\n#####################\n\n')
fprintf(['The resulting residuals for both the polynomial and linear functions are presented down below:\n\n' ...
    'Epsilon_polynomial = %.f\nEpsilon_linear = %.f\n\n'], eps_curve, eps_linear);



%% 2.3a)
% Load and defining variables
u = load("input.mat");
y = load("output.mat");
N = length(u.u);
var = 0.01;
e =  normrnd(0, var, [N 1]);

%% 3-(12a)
% Dividing set into pred and val set
uest = u.u(1:N/2);
yest = y.y(1:N/2);
uval = u.u(N/2+1:end);
yval = y.y(N/2+1:end);

% H for 12a)
H = zeros(N/2,3);
H(1,:) = [0       0    uest(1)];
H(2,:) = [-yest(1) 0    uest(2)];
for t=3:N/2
    H(t,:) = [-yest(t-1) -yest(t-2) uest(t)];
end

% Calculating belonging parameters
th_hat_a = (H'*H)\H'*yest;
a0_hat_a = th_hat_a(1);
a1_hat_a = th_hat_a(2);
b0_hat_a = th_hat_a(3);

% Prediction step for 12a)
NN = N/2;
yn = yval;
un = uval;
ypred = zeros(NN,1); 
ypred(1) = yn(1);
ypred(2) = yn(2);
for t=3:NN
     ypred(t) = [-yn(t-1) -yn(t-2) un(t)]*[a0_hat_a; a1_hat_a; b0_hat_a];
end

% Calculating and printing RMSE for 12a) pred
predERROR = yn-ypred;
predRMSE  = rms(predERROR);
fprintf('\n#####################\n     Model 12a)\n#####################\n')
disp(['a1: ', num2str(a0_hat_a), ', a2: ', num2str(a1_hat_a) ...
    ', b0: ', num2str(b0_hat_a)])
disp(['Prediction RMS error for (12a) is: ' num2str(predRMSE)])

% Simulation step for 12a)
ysim = zeros(NN,1); 
ysim(1) = yn(1); 
ysim(2) = yn(2); 
for t=3:NN 
    ysim(t) = -a0_hat_a * ysim(t-1) - a1_hat_a * ysim(t-2)  + b0_hat_a*un(t-1);
end

% Calculating and printing RMSE for 12a) sim
simERROR = yn-ysim;
simRMSE = rms(simERROR);
disp(['Simulation RMS error for (12a) is: ' num2str(simRMSE)])



%% 3-(12b)
% H for 12b)
H = zeros(N/2,4);
H(1,:) = [-yest(1) yest(1) uest(1) uest(1)];
H(2,:) = [-yest(1) yest(1) uest(2) uest(1)];
for t=3:N/2
    H(t,:) = [-yest(t-1) -yest(t-2) uest(t) uest(t-1)];
end

% Calculating belonging parameters
th_hat_b = (H'*H)\H'*yest;
a0_hat_b = th_hat_b(1);
a1_hat_b = th_hat_b(2);
b0_hat_b = th_hat_b(3);
b1_hat_b = th_hat_b(4);

% Prediction step for 12b)
ypred = zeros(NN,1); 
ypred(1) = yn(1);
ypred(2) = yn(2);
for t=3:NN 
     ypred(t) = [-yn(t-1) -yn(t-2) un(t) un(t-1)]*[a0_hat_b; a1_hat_b; b0_hat_b; b1_hat_b];
end

% Calculating and printing RMSE for 12b) pred
predERROR = yn-ypred;
predRMSE  = rms(predERROR);
fprintf('\n#####################\n     Model 12b)\n#####################\n')
disp(['a1: ', num2str(a0_hat_b), ', a2: ', num2str(a1_hat_b) ...
    ' b0: ', num2str(b0_hat_b) ', b1: ', num2str(b1_hat_b)])
disp(['Prediction RMS error for (12b) is: ' num2str(predRMSE)])

% Simulation step for 12b)
ysim = zeros(NN,1); 
ysim(1) = yn(1); 
ysim(2) = yn(2); 
for t=3:NN 
    ysim(t) = -a0_hat_b * ysim(t-1) - a1_hat_b * ysim(t-2)  + b0_hat_b*un(t) + b1_hat_b*un(t-1);
end

% Calculating and printing RMSE for 12b) sim
simERROR = yn-ysim;
simRMSE = rms(simERROR);
disp(['Simulation RMS error for (12b) is: ' num2str(simRMSE)])



%% 3-(12c)
% H for 12c)
H = zeros(N/2,4);
H(1,:) = [ 0           0        0    0];
H(2,:) = [-yest(1)     0        0    uest(1)];
H(3,:) = [-yest(2)  -yest(1)     0   uest(2)];
for t=4:N/2
    H(t,:) = [-yest(t-1) -yest(t-2) -yest(t-3) uest(t-1)];
end

% Calculating belonging parameters
th_hat_c = (H'*H)\H'*yest;
a0_hat_c = th_hat_c(1);
a1_hat_c = th_hat_c(2);
a2_hat_c = th_hat_c(3);
b1_hat_c = th_hat_c(4);

% Prediction step for 12c)
ypred = zeros(NN,1); 
ypred(1) = yn(1);
ypred(2) = yn(2);
ypred(3) = yn(3);
for t=4:NN
     ypred(t) = [-yn(t-1) -yn(t-2) -yn(t-3) un(t-1)]*[a0_hat_c; a1_hat_c; a2_hat_c; b1_hat_c];
end

% Calculating and printing RMSE for 12c) pred
predERROR = yn-ypred;
predRMSE  = rms(predERROR);
fprintf('\n#####################\n     Model 12c)\n#####################\n')
disp(['a1: ', num2str(a0_hat_c), ', a2:', num2str(a1_hat_c) ...
    ', a3: ', num2str(a2_hat_c) ', b1: ', num2str(b1_hat_c)])
disp(['Prediction RMS error for (12c) is: ' num2str(predRMSE)])

% Simulation step for 12b)
ysim = zeros(NN,1); 
ysim(1) = yn(1); 
ysim(2) = yn(2); 
ysim(3) = yn(3);
for t=4:NN 
    ysim(t) = -a0_hat_c * ysim(t-1) - a1_hat_c * ysim(t-2) - a2_hat_c*ysim(t-3) + b1_hat_c*un(t-1);
end

% Calculating and printing RMSE for 12c) sim
simERROR = yn-ysim;
simRMSE = rms(simERROR); 
disp(['Simulation RMS error for (12c) is: ' num2str(simRMSE)])


