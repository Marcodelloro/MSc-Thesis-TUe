% Define the symbolic variables
clc
clear all
close all

syms S I D T1 T2 alpha beta gamma delta1 delta2 sigma1 tau1 sigma2 tau2 epsilon lambda
syms S I D T H alpha gamma delta sigma tau lambda

% % Define the Jacobian matrix J_F_X
% J_F_X = [S*alpha, S*beta, 0, 0;
%          0,       0,      0, 0;
%          0,       0,      0, 0;
%          0,       0,      0, 0];
% 
% % Define the Jacobian matrix J_V_X
% J_V_X = [gamma + lambda,                   0, 0,                     0;
%          -gamma,          delta1 + delta2 + lambda, 0,                     0;
%          0,                -delta1,          epsilon + sigma1 + tau1, 0;
%          0,                -delta2,          epsilon,                 sigma2 + tau2];

% Define the Jacobian matrix J_F_X
J_F_X = [   S*alpha, 0,    0;
            0,       0,    0;
            0,       0,    0;
            0,       0,    0  ];

% Define the Jacobian matrix J_V_X
J_V_X = [   gamma + lambda,        0,               0;
             -gamma,          delta + lambda,       0;
                0,              -delta,         sigma + tau2     ];


r0 = J_F_X * inv(J_V_X)

%% Include the symbolic math toolbox

syms S S_dot I I_dot D D_dot T1 T1_dot T2 T2_dot H H_dot E V alpha beta gamma lambda delta1 delta2 sigma1 tau1 epsilon sigma2 tau2 k V_psi psi real

% Define the equations

Eqs = [ -S*(alpha + beta*D) - S*psi + V*psi;
         S*(alpha + beta*D) - (gamma + lambda)*I;
         I*gamma - D*(lambda + delta1 + delta2);
         delta1*D - (sigma1 + tau1 + epsilon)*T1;
         delta2*D - (sigma2 + tau2)*T2 + epsilon*T1;
         (I + D)*lambda + T1*sigma1 + T2*sigma2 - H*k;
         T1*tau1 + T2*tau2;
         S*psi + H*k - V*psi ];

x = [S, I, D, T1, T2, H, E, V];

A = jacobian(Eqs(:), x)

syms c a beta

A = [-c -a;
        1 0]
eig(A)

%%
syms S I D T H alpha gamma delta sigma tau lambda p

J = [p,             -S*alpha,             0,                0,           0;
     0,     p - S*alpha + gamma + lambda, 0,                0,           0;
     0,               gamma,         p + delta + lambda,    0,           0;
     0,                 0,               delta,      p + tau + sigma,  0;
     0,             lambda,             lambda,          sigma,          p      ];


a = det(J);
eq = collect(a,p)

solve(eq == 0, p)

%% Camputation of Jacobian od SIDTTHE

syms S I D T1 T2 H E V alpha beta gamma lambda delta1 delta2 sigma1 sigma2 tau1 tau2 epsilon phi psi kappa p

% Define the equations
S_dot = -S*(alpha*I + beta*D) - S*phi + V*psi;
I_dot = S*(alpha*I + beta*D) - (gamma + lambda)*I;
D_dot = I*gamma - D*(lambda + delta1 + delta2);
T1_dot = delta1*D - (sigma1 + tau1 + epsilon)*T1;
T2_dot = delta2*D - (sigma2 + tau2)*T2 + epsilon*T1;
H_dot = (I + D)*lambda + T1*sigma1 + T2*sigma2-kappa*H;
E_dot = T1*tau1 + T2*tau2;
V_dot = S*phi + kappa*H - V*psi;

% Create a vector of these equations
equations = [S_dot; I_dot; D_dot; T1_dot; T2_dot; H_dot; E_dot; V_dot];

% Calculate the Jacobian matrix of the system with respect to the variables
variables = [S I D T1 T2 H E V];
jacobian_matrix = jacobian(equations, variables);
jac_det = jacobian_matrix*-1

J2 = [    p+phi,                  S*alpha,                       S*beta,                       0,             0,            0,    0,        -psi;
            0,           p + gamma + lambda - S*alpha,          -S*beta,                       0,             0,            0,    0,          0;
            0,                   -gamma,              p + delta1 + delta2 + lambda,            0,             0,            0,    0,          0;
            0,                        0,                        -delta1,        p + epsilon + sigma1 + tau1,  0,            0,    0,          0;
            0,                        0,                        -delta2,                   -epsilon,   p + sigma2 + tau2,   0,    0,          0;
            0,                    -lambda,                      -lambda,                    -sigma1,      -sigma2,     p + kappa, 0,          0;
            0,                        0,                           0,                        -tau1,         -tau2,          0,    p,          0;
          -phi,                       0,                           0,                          0,             0,         -kappa,  0,         p+psi   ];

% J2 = [    p+phi,                  S*alpha,                       S*beta,                       0,             0,            0,       0;
%             0,           p + gamma + lambda - S*alpha,          -S*beta,                       0,             0,            0,       0;
%             0,                    -gamma,             p + delta1 + delta2 + lambda,            0,             0,            0,       0;
%             0,                        0,                        -delta1,        p + epsilon + sigma1 + tau1,  0,            0,       0;
%             0,                        0,                        -delta2,                   -epsilon,   p + sigma2 + tau2,   0,       0;
%             0,                    -lambda,                      -lambda,                    -sigma1,      -sigma2,        p+kappa,   0;
%            -phi,                       0,                           0,                          0,             0,         -kappa,    p  ];

a = det(J2);
eq = collect(a,p)

solve(eq == 0, p)

%% Camputation of Jacobian od SIDTHE

syms S I D T H E alpha gamma lambda delta sigma tau p

% Define the equations
S_dot = -S*(alpha*I);
I_dot = S*(alpha*I) - (gamma + lambda)*I;
D_dot = I*gamma - D*(lambda + delta);
T_dot = delta*D - (sigma + tau)*T;
H_dot = (I + D)*lambda + T*sigma;
E_dot = T*tau;

% Create a vector of these equations
equations = [S_dot; I_dot; D_dot; T_dot; H_dot; E_dot];

% Calculate the Jacobian matrix of the system with respect to the variables
variables = [S I D T H E];
jacobian_matrix = jacobian(equations, variables)
jac_det = jacobian_matrix*-1

J2 =        [ p + I*alpha,                  S*alpha,                 0,             0,          0,          0;
                -I*alpha,         p + gamma + lambda - S*alpha,      0,             0,          0,          0;
                    0,                      -gamma,        p + delta + lambda,      0,          0,          0;
                    0,                          0,                -delta,   p + sigma + tau,    0,          0;
                    0,                      -lambda,              -lambda,       -sigma,        p,          0;
                    0,                          0,                   0,           -tau,         0,          p   ];



a = det(J2);
eq = collect(a,p)

solve(eq == 0, p)

