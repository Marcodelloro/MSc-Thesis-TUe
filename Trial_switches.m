% Trial To make a set of two equations switch based on a 

clc
clear all

addpath('/Users/marcodelloro/Downloads/casadi-3.6.3-osx64-matlab2018b')
import casadi.*;
opti = casadi.Opti();

%% Initialisation

% Collocation settings
problem_options = NosnocProblemOptions();
solver_options = NosnocSolverOptions();
model = NosnocModel();

problem_options.irk_scheme = IRKSchemes.RADAU_IIA;
problem_options.cross_comp_mode = 1;
problem_options.n_s = 2;

% Time-settings  - Solve an time optimal control problem
problem_options.time_optimal_problem = 1;

% Penalty/Relaxation paraemetr
solver_options.comp_tol = 1e-9;
problem_options.cross_comp_mode = 1;
problem_options.N_stages = 10;         % Number of control intervals
problem_options.N_finite_elements = 3; % Number of finite element on every control interval (optionally a vector might be passed)
problem_options.T = 1;                 % Time horizon

%% Info

% The goal of to make the system of two equation switch as they cross a discontinuty (the value of T2 becomes too big)
% Here is created a model of the SIDTTHE with just T1 and T2

N = 1;

delta1 = 0.01 ;   % delta 1 - coefficient to go from D to T1
delta2 = 0.001 ;  % delta 2 - coefficient to go from D to T2 
epsi = 0.002 ;    % epsilon - coefficient to go from T1 to T2
sigma1 = 0.002 ;  % sigma 1 - coefficient to go from T1 to H

% delta1 = 0.01 * ones(1,N);  % delta 1 - coefficient to go from D to T1
% delta2 = 0.001 * ones(1,N);  % delta 2 - coefficient to go from D to T2 
% epsi = 0.002 * ones(1,N);   % epsilon - coefficient to go from T1 to T2
% sigma1 = 0.002 * ones(1,N);  % sigma 1 - coefficient to go from T1 to H

coefs = [delta1; delta2; epsi; sigma1];

% Lifting Variables

T1 = opti.variable(1,N);
T2 = opti.variable(1,N);

X = [T1; T2];      % state trajectory

opti.subject_to(T1 >= 0);
opti.subject_to(T2 >= 0);

ICUs = 2000;

%% Specify system dynamics - creation of the handle functions

% f = @(X, coefs)   [          ( (coefs(1) + coefs(2) + coefs(3))*X(1) );                            % dT1/dt = T1_dot
%                              ( (coefs(4) + coefs(3))*X(1) ) + ( coefs(1)*X(2) )                    % dT2/dt = T2_dot
%                   ];
% 
% dt = 1;
% for k=1:N-1      % loop over control intervals
% 
% % Runge-Kutta 4 integration
% 
%    k1 = f(X(:,k),         coefs(:,k));
%    k2 = f(X(:,k)+dt/2*k1, coefs(:,k));
%    k3 = f(X(:,k)+dt/2*k2, coefs(:,k));
%    k4 = f(X(:,k)+dt*k3,   coefs(:,k));
%    x_next = X(:,k) + dt/6*(k1+2*k2+2*k3+k4);
%    opti.subject_to(X(:,k+1) == x_next); % close the gaps
% 
% end

%% System Switch

% Solver options

solver_options = NosnocSolverOptions();
solver_options.store_integrator_step_results = 1;
problem_options.T_sim = 1.5;
problem_options.N_sim = 7;

model.x = X ;
model.c = (coefs(4) + coefs(3))*X(1) + ( coefs(1)*X(2) ) - ICUs ;
model.S = [1; -1] ;

f_1 = [      ( (coefs(1) + coefs(2) + coefs(3))*X(1) );                            % dT1/dt = T1_dot
             ( (coefs(4) + coefs(3))*X(1) ) + ( coefs(1)*X(2) )                    % dT2/dt = T2_dot
      ];

f_2 = [       ( (coefs(1)*coefs(3) - coefs(2) + coefs(3))*X(1) );                  % dT1/dt = T1_dot
              ( (coefs(4) + coefs(3))*X(1) ) + ( coefs(1)*X(2)*X(1) )              % dT2/dt = T2_dot
      ];

model.F = [f_1 f_2] ;      % Matrix which states the modes of f(x)

%% Solve OCP

mpcc = NosnocMPCC(problem_options, model);
solver = NosnocSolver(mpcc, solver_options);
[results,stats] = solver.solve();

