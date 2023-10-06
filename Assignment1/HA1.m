clear all; close all; clc

% Define symbolic variables
syms p1 theta phi real;
syms  real; % Define p2 as a function of q and f
syms m2 g real; % Mass of the hanging object and gravitational constant
syms u1 u2 u3 real; % External force applied to the helicopter

% Define generalized coordinates
q = [p1; theta; phi];

% Calculate position vector of p2
p2 = p1 + [l*sin(theta)*cos(phi); l*sin(theta)*sin(phi); -l*cos(theta)];

% Calculate derivatives of p2 with respect to q and f
dp2_dq = diff(p2, theta);
dp2_df = diff(p2, phi);

% Calculate kinetic energy T2 and potential energy V2
T2 = 0.5 * m2 * (dp2_dq' * dp2_dq + dp2_df' * dp2_df);
V2 = m2 * g * [0; 0; 1] * p2;

% Calculate Lagrangian L
L = T2 - V2;

% Compute the Lagrange equations
% Calculate the generalized forces
Q = [u1; u2; u3]; % Assuming u is applied to p1

% Calculate the equations of motion (Lagrange's equations)
syms p1_dot q_dot f_dot real; % Define generalized velocities
q_dot = [p1_dot; q_dot; f_dot];
eqns = simplify(jacobian(jacobian(L, q_dot), q) - jacobian(L, q)' - Q);

% Define the mass matrix M and generalized forces vector b
M = jacobian(eqns, q_dot);
b = simplify(eqns - M * q_dot);

% Now you have M(q) and b(q, q_dot, u), which represent your dynamic model

% Define initial conditions and external forces
q0 = [0; 0; 0]; % Initial generalized coordinates
q_dot0 = [0; 0; 0]; % Initial generalized velocities
u = [u1; u2; u3]; % External forces applied to the helicopter

% Solve for the accelerations q_ddot
q_ddot = simplify(M \ (b + M * q_dot0));

% Example: Evaluate q_ddot for specific conditions
q_ddot_values = subs(q_ddot, {q, q_dot, u}, {q0, q_dot0, u});

% Display the result
disp('q_ddot for given conditions:');
disp(q_ddot_values);
