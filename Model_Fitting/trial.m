%% Version of V4 - V3 but coded better

clc
clear all
close all

load('/Users/marcodelloro/Desktop/Thesis/MSc-Thesis-TUe/Data_Collection/SIDTTHE_data.mat');
load('/Users/marcodelloro/Desktop/Thesis/MSc-Thesis-TUe/Data_Collection/Tests_data.mat');
load('/Users/marcodelloro/Desktop/Thesis/MSc-Thesis-TUe/Data_Collection/vaxData.mat'); 

coef_val = load('/Users/marcodelloro/Desktop/matfile.mat'); 

addpath('/Users/marcodelloro/Downloads/casadi-3.6.3-osx64-matlab2018b')
import casadi.*;
opti = casadi.Opti();

[DataArray, new_time, matrixData] = LoadingData(SIDTTHE_data,vaxData);
save("MartixDataSIDDTHE.mat", 'matrixData');
data = cell2mat(DataArray); 
N = 399; % number of days

%% Setting of the optimization variables

alpha = opti.variable(1,1);   % alpha - coefficient to go from S to I 
beta = opti.variable(1,1);    % beta - coefficient to go from S to D
gamma = opti.variable(1,1);   % gamma - coefficient to go from I to D
delta1 = opti.variable(1,1);  % delta 1 - coefficient to go from D to T1
delta2 = opti.variable(1,1);  % delta 2 - coefficient to go from D to T2 
epsi = opti.variable(1,1);    % epsilon - coefficient to go from T1 to T2
sigma1 = opti.variable(1,1);  % sigma 1 - coefficient to go from T1 to H
sigma2 = opti.variable(1,1);  % sigma 2 - coefficient to go from T2 to H
tau1 = opti.variable(1,1);    % tau 1 - coefficient to go from T1 to E
tau2 = opti.variable(1,1);    % tau 2 - coefficient to go from T2 to E
lambda = opti.variable(1,1);  % lambda - coefficient to go from I and D to H
phi = opti.variable(1,1);     % phi - coefficient to go from S to V
kappa = opti.variable(1,1);   % phi - coefficient to go from H to V

coefs = [alpha; beta; gamma; delta1; delta2; epsi; sigma1; sigma2; tau1; tau2; lambda; phi; kappa];

% Lifting Variables

S = opti.variable(1,1); 
I = opti.variable(1,1);
D = opti.variable(1,1);
T1 = opti.variable(1,1);
T2 = opti.variable(1,1);
H = opti.variable(1,1);
E = opti.variable(1,1);
V = opti.variable(1,1);

X = [S; I; D; T1; T2; H; E; V]; % state trajectory

% Initial conditions for the various parameters

% opti.set_initial(S,ones(N,1) * data(1,1));
% opti.set_initial(I,ones(N,1) * data(2,1));
% opti.set_initial(D,ones(N,1) * data(3,1));
% opti.set_initial(T1,ones(N,1)* data(4,1));
% opti.set_initial(T2,ones(N,1)* data(5,1));
% opti.set_initial(H,ones(N,1) * data(6,1));
% opti.set_initial(E,ones(N,1) * data(7,1));
% opti.set_initial(V,ones(N,1) * data(8,1));
% 
% opti.set_initial(alpha,ones(N,1) * 0.513);
% opti.set_initial(beta,ones(N,1) * 0.011);
% opti.set_initial(beta,ones(N,1) * 1e-6);
% opti.set_initial(gamma,ones(N,1) * 0.2405);
% opti.set_initial(delta1,ones(N,1) * 0.005);
% opti.set_initial(delta2,ones(N,1) * 0.001);
% 
% opti.set_initial(sigma1,ones(N,1) * 0.03);
% opti.set_initial(sigma2,ones(N,1) * 0.065);
% opti.set_initial(tau1,ones(N,1) * 0.015);
% opti.set_initial(tau2,ones(N,1) * 0.005);
% opti.set_initial(lambda,ones(N,1) * 0.0596);

%% Constraints for the parameters
% This part is the one that presents most of the difference with respect to
% the previous model; as a matter of facts, we do not have a cost to
% smoothing the parameters, but constraints

% Alpha parameter - bound on alfa value (different bound for different months)

opti.subject_to(alpha >= 0);
opti.subject_to(alpha <= 0.6); 

% Beta parameter - bound on beta value (different bound for different months)

opti.subject_to(beta(1:end) >= 0);
opti.subject_to(beta(1:end) <= 0.015);

% Gamma parameter - bound on gamma value (different bound for different months)

opti.subject_to(gamma >= 0); 
opti.subject_to(gamma <= 0.5);

% delta parameter - bound on delta value (different bound for different months)

opti.subject_to(delta1 >= 0);     % bound on delta1 value
opti.subject_to(delta1 <= 0.002);

opti.subject_to(delta2 >= 0);     % bound on delta2 value
opti.subject_to(delta2 <= 5e-4);

% epsi parameter - bound on epsilon value (different bound for different months)

opti.subject_to(epsi >= 0);    % bound on epsi value
opti.subject_to(epsi <= 5e-4);

% sigma parameter - bound on sigma value (different bound for different months)

opti.subject_to(sigma1 >= 0);   % bound on sigma1 value
opti.subject_to(sigma1 <= 0.03);

opti.subject_to(sigma2 >= 0);   % bound on sigma2 value
opti.subject_to(sigma2 <= 0.03);

% tau parameter - bound on tau value (different bound for different months)

opti.subject_to(tau1 >= 0);     % bound on tau1 value
opti.subject_to(tau1 <= 0.03);

opti.subject_to(tau2 >= 0);     % bound on tau2 value
opti.subject_to(tau2 <= 0.02);

% Lambda parameter - bound on lambda value (different bound for different months)

opti.subject_to(lambda >= 0);   % bound on lambda value
opti.subject_to(lambda <= 0.2);

% phi parameter - bound on phi value (different bound for different months)

opti.subject_to(phi >= 0);     % bound on phi value
opti.subject_to(kappa >= 0);   % bound on phi value

% Initial Conditions constraints on ODEs

% % opti.subject_to(S(1,1) == data(1,1));
% % opti.subject_to(I(1,1) == data(2,1));
% % opti.subject_to(D(1,1) == data(3,1));
% % opti.subject_to(T1(1,1) == data(4,1));
% % opti.subject_to(T2(1,1) == data(5,1));
% % opti.subject_to(H(1,1) == data(6,1));
% % opti.subject_to(E(1,1) == data(7,1));
% % opti.subject_to(V(1,1) == data(8,1));
% 
% % Initial Conditions constraints on Parameters 
% 
% % opti.subject_to(alpha(1,1) == 0.2);
% % opti.subject_to(beta(1,1) == 10e-3);
% % opti.subject_to(gamma(1,1) == 0.3);
% 
opti.subject_to(S + I + D + T1 + T2 + H + E + V == 1);
% 
% opti.subject_to(S >= 0);
% opti.subject_to(I >= 0);
% opti.subject_to(D >= 0);
% opti.subject_to(T1 >= 0);
% opti.subject_to(T2 >= 0);
% opti.subject_to(H >= 0);
% opti.subject_to(E >= 0);
% opti.subject_to(V >= 0);
% 
% opti.subject_to(alpha >= 2*beta); % bound on contagion parameter between alpha and beta
% 
% % my constraints assumptions 
% 
% opti.subject_to(tau1 <= tau2);
% opti.subject_to(sigma2 <= sigma1);
% opti.subject_to(delta1 > delta2);

%% HARD Constraints on the variation of the coefficents 
% 
% % Constraints related to the coefficients - Note that coefficients smoothness will be required only in
% % 'constant zones' so in time perioods with constant policies
% % For subsequent "dd" value, the constraints will be removed 
% 
% 8 different areas will be covered
policy_idx = [1 40 69 116 141 190 242 294 399];   % Values taken from italian policies applied for the pandemic

for ii=1:length(policy_idx)
    policy_dates(ii) = [new_time(policy_idx(ii))];
end
% 
% dd = 1; % days for how much the constraint will be removed

% Note that as assumption I ASSUME that only coefficients \alpha and \beta can be influenced by
% the different policies applied form the government

% constraints NO VARIATION 
% for kk = 1:length(policy_idx)-1
%     for jj = policy_idx(kk)+dd:policy_idx(kk+1)-1
% 
%         opti.subject_to( alpha(1,jj + 1) == alpha(1,jj) )
%         opti.subject_to( alpha(1,jj + 1) == alpha(1,jj) )
% 
%         opti.subject_to( beta(1,jj + 1) == beta(1,jj) )
%         opti.subject_to( beta(1,jj + 1) == beta(1,jj) )
% 
%     end
% end

%---------------------------------------------------------------------

% % constraints WITH VARIATION 
% for kk = 2:length(policy_idx)-1
%     for jj = policy_idx(kk)+dd:policy_idx(kk+1)-1
% 
%     opti.subject_to( alpha(1,jj + 1) <= alpha(1,jj)*1.01 )
%     opti.subject_to( alpha(1,jj + 1) >= alpha(1,jj)*0.99 )
% 
%     opti.subject_to( beta(1,jj + 1) <= beta(1,jj)*1.01 )
%     opti.subject_to( beta(1,jj + 1) >= beta(1,jj)*0.99 )
% 
%     end
% end

%% Additional constraints to make the Coefs Work

% constraint on the first part of the simulation
% for jj = 1:policy_idx(2)-1
%     opti.subject_to( alpha(1,jj + 1) <= alpha(1,jj)*1.01 )
%     opti.subject_to( alpha(1,jj + 1) >= alpha(1,jj)*0.99 )
% 
%     % opti.subject_to( beta(1,jj + 1) <= beta(1,jj)*1.01 )
%     % opti.subject_to( beta(1,jj + 1) >= beta(1,jj)*0.99 )
% end

%% RK45 simulation of the differential equations

% Since the equations are always the same, we can just load the already solved RK45

rhs =               [ -X(1)*(coefs(1)*X(2) + coefs(2)*X(3)) - X(1)*coefs(12) ;                              % dS/dt = S_dot 
                    ( X(1)*(coefs(1)*X(2) + coefs(2)*X(3)) ) - ( (coefs(3) + coefs(11))*X(2) );             % dI/dt = I_dot
                    ( X(2)*coefs(3)) - ( X(3)*(coefs(11) + coefs(4) + coefs(5)) );                          % dD/dt = D_dot
                    ( coefs(4)*X(3) ) - ( (coefs(7) + coefs(9) + coefs(6))*X(4) );                          % dT1/dt = T1_dot
                    ( coefs(5)*X(3) ) - ( (coefs(8) + coefs(10))*X(5) ) + ( coefs(6)*X(4) );                % dT2/dt = T2_dot
                    ( (X(2) + X(3))*coefs(11) ) + ( X(4)*coefs(7) ) + ( X(5)*coefs(8) ) - X(6)*coefs(13) ;  % dH/dt = H_dot
                    ( X(4)*coefs(9) ) + ( X(5)*coefs(10) );
                     X(1)*coefs(12) + X(6)*coefs(13) ];
     
ode = Function('ode',{X,coefs},{rhs});

dt = 1;

% Build an integrator for this system: Runge Kutta 4 integrator

k1 = ode(X,coefs);
k2 = ode(X+dt/2.0*k1,coefs);
k3 = ode(X+dt/2.0*k2,coefs);
k4 = ode(X+dt*k3,coefs);
states_final = X+dt/6.0*(k1+2*k2+2*k3+k4);

% Create a function that simulates one step propagation in a sample
one_step = Function('one_step',{X,coefs},{states_final});

X_sim = X;
for i=1:1
    X_sim = one_step(X_sim, coefs);
end

% Create a function that simulates all step propagation on a sample
one_sample = Function('one_sample',{X, coefs}, {X_sim});
one_sample = one_sample.expand();

x0 = DM([data(:,1)]);

% Simulating the system 
all_samples = one_sample.mapaccum('all_samples', N);
X_symbolic = all_samples(x0, repmat(coefs,1,N));
                             
%% Actual Optimization Simulation

data_obj = horzcat(data(:));
X_obj = horzcat(X_symbolic(:));
weightVec = [0.5 0.25 0.75 1 1 0.25 0.1 0.75];
Qvec = repmat(weightVec',399,1)';

maxData =( ones(399,1)*max(data,[],2)' )';
maxData = horzcat(maxData(:));

meanData =( ones(399,1)*mean(data,2)' )';
meanData = horzcat(meanData(:));

% normalization in cost function - max data

obj = sum(((data_obj - X_obj)./maxData).^2);

opti.minimize(obj);
p_opts = struct('expand', false);
s_opts = struct('max_iter',1e4, 'tol',1e-4,'constr_viol_tol',1e-4,'compl_inf_tol',1e-4,'linear_solver','MA97'); %'MA27' 'MA57' 'MA77' 'MA86' 'MA97'
opti.solver('ipopt',p_opts, s_opts); % set numerical backend
sol = opti.solve();   % actual solver

%% Gathering of all data

% Table of the coefficients
column_names = {'S', 'I', 'D', 'T1', 'T2', 'H', 'E', 'V'};
SIDTTHE_trends = array2table(reshape(opti.debug.value(X_obj),[8,399])', 'VariableNames', column_names);


column_names = {'alpha', 'beta', 'gamma', 'delta1', 'delta2', 'epsi', 'sigma1', 'sigma2', 'tau1', 'tau2', 'lambda','phi','kappa'};
opti_coefficients = array2table(opti.debug.value(coefs)', 'VariableNames', column_names);

% Trends of every group SIDTTHE
% column_names = {'S', 'I', 'D', 'T1', 'T2', 'H', 'E', 'V'};
% SIDTTHE_trends = array2table(opti.debug.value(X)', 'VariableNames', column_names);

% % check if the coefficients reach the bounds
% 
% alpha_bound = any(opti_coefficients.alpha == 0.6)
% beta_bound = any(opti_coefficients.beta == 0.015)
% gamma_bound = any(opti_coefficients.gamma == 0.5)
% delta1_bound = any(opti_coefficients.delta1 == 0.002)
% delta2_bound = any(opti_coefficients.delta2 == 5e-4)
% tau1_bound = any(opti_coefficients.tau1 == 0.03)
% tau2_bound = any(opti_coefficients.tau2 == 0.05)
% sigma1_bound = any(opti_coefficients.sigma1 == 0.03)
% sigma2_bound = any(opti_coefficients.sigma2 == 0.03)
% lambda_bound = any(opti_coefficients.lambda == 0.2)
% 
% % Saving of the simulation results 
% 
% Opti_results = {SIDTTHE_trends, 'Group Trends'; opti_coefficients, 'Coefficients Optimised'; new_time, 'Time for Plotting'};
% filename = 'Optimization_Results.mat';
% save(filename, 'Opti_results');

%% Comparison Plot - ODE curves

print_curves(SIDTTHE_trends,DataArray,new_time)

% %% Comparison Plot - Coefficients
% load("/Users/marcodelloro/Desktop/Thesis/MSc-Thesis-TUe/Data_Collection/Smooth_Variants.mat")
% 
% print_coefs(opti_coefficients,new_time, policy_dates,SmoothVar,tests)
