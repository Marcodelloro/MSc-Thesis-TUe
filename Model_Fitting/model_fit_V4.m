%% Version of V4 - V3 but coded better

clc
clear all
close all

load('/Users/marcodelloro/Desktop/Thesis/MSc-Thesis-TUe/Data_Collection/SIDTTHE_data.mat');
load('/Users/marcodelloro/Desktop/Thesis/MSc-Thesis-TUe/Data_Collection/Tests_data.mat');
load('/Users/marcodelloro/Desktop/Thesis/MSc-Thesis-TUe/Data_Collection/vaxData.mat'); 
addpath('/Users/marcodelloro/Downloads/casadi-3.6.3-osx64-matlab2018b')
import casadi.*;
opti = casadi.Opti();

[DataArray, new_time, matrixData] = LoadingData(SIDTTHE_data,vaxData);
save("MartixDataSIDDTHE.mat", 'matrixData');
data = cell2mat(DataArray); 

% Estimation of parameters from data up to a certain value - 30 november 2020
N = 399;
shortData = data(1:7,1:N);


%% Setting of the optimization variables

alpha = opti.variable(1,N);   % alpha - coefficient to go from S to I 
beta = opti.variable(1,N);    % beta - coefficient to go from S to D
gamma = opti.variable(1,N);   % gamma - coefficient to go from I to D
delta1 = opti.variable(1,N);  % delta 1 - coefficient to go from D to T1
delta2 = opti.variable(1,N);  % delta 2 - coefficient to go from D to T2 
epsi = opti.variable(1,N);   % epsilon - coefficient to go from T1 to T2
sigma1 = opti.variable(1,N);  % sigma 1 - coefficient to go from T1 to H
sigma2 = opti.variable(1,N);  % sigma 2 - coefficient to go from T2 to H
tau1 = opti.variable(1,N);    % tau 1 - coefficient to go from T1 to E
tau2 = opti.variable(1,N);    % tau 2 - coefficient to go from T2 to E
lambda = opti.variable(1,N);  % lambda - coefficient to go from I and D to H

coefs = [alpha; beta; gamma; delta1; delta2; epsi; sigma1; sigma2; tau1; tau2; lambda];

% Lifting Variables

S = opti.variable(1,N); 
I = opti.variable(1,N);
D = opti.variable(1,N);
T1 = opti.variable(1,N);
T2 = opti.variable(1,N);
H = opti.variable(1,N);
E = opti.variable(1,N);

X = [S; I; D; T1; T2; H; E]; % state trajectory

% Initial conditions for the various parameters

opti.set_initial(S,ones(N,1) * data(1,1));
opti.set_initial(I,ones(N,1) * data(2,1));
opti.set_initial(D,ones(N,1) * data(3,1));
opti.set_initial(T1,ones(N,1)* data(4,1));
opti.set_initial(T2,ones(N,1)* data(5,1));
opti.set_initial(H,ones(N,1) * data(6,1));
opti.set_initial(E,ones(N,1) * data(7,1));

opti.set_initial(alpha,ones(N,1) * 0.513);          % All coeff. are obtained from the mean of the values provided by G.Giordano in original SIDARTHE paper
opti.set_initial(beta,ones(N,1) * 0.011);
opti.set_initial(beta,ones(N,1) * 1e-6);
opti.set_initial(gamma,ones(N,1) * 0.2405);
opti.set_initial(delta1,ones(N,1) * 0.005);
opti.set_initial(delta2,ones(N,1) * 0.001);

opti.set_initial(sigma1,ones(N,1) * 0.03);
opti.set_initial(sigma2,ones(N,1) * 0.065);
opti.set_initial(tau1,ones(N,1) * 0.015);
opti.set_initial(tau2,ones(N,1) * 0.005);
opti.set_initial(lambda,ones(N,1) * 0.0596);

%% Constraints for the parameters
% This part is the one that presents most of the difference with respect to
% the previous model; as a matter of facts, we do not have a cost to
% smoothing the parameters, but constraints

% Alpha parameter - bound on alfa value (different bound for different months)

opti.subject_to(alpha(1:end) >= 0);
opti.subject_to(alpha(1:end) <= 0.6); 

% Beta parameter - bound on beta value (different bound for different months)

opti.subject_to(beta(1:end) >= 0);
opti.subject_to(beta(1:end) <= 0.015);

% Gamma parameter - bound on gamma value (different bound for different months)

opti.subject_to(gamma(1:end) >= 0); 
opti.subject_to(gamma(1:end) <= 0.5);

% delta parameter - bound on delta value (different bound for different months)

opti.subject_to(delta1(1:end) >= 0);     % bound on delta1 value
opti.subject_to(delta1(1:end) <= 0.002);

opti.subject_to(delta2(1:end) >= 0);     % bound on delta2 value
opti.subject_to(delta2(1:end) <= 5e-4);

% epsi parameter - bound on epsilon value (different bound for different months)

opti.subject_to(epsi(1:end) >= 0);    % bound on epsi value
opti.subject_to(epsi(1:end) <= 5e-4);

% sigma parameter - bound on sigma value (different bound for different months)

opti.subject_to(sigma1(1:end) >= 0);   % bound on sigma1 value
opti.subject_to(sigma1(1:end) <= 0.03);

opti.subject_to(sigma2(1:end) >= 0);   % bound on sigma2 value
opti.subject_to(sigma2(1:end) <= 0.03);

% tau parameter - bound on tau value (different bound for different months)

opti.subject_to(tau1(1:end) >= 0);     % bound on tau1 value
opti.subject_to(tau1(1:end) <= 0.03);

opti.subject_to(tau2(1:end) >= 0);     % bound on tau2 value
% opti.subject_to(tau2(1:54) <= 0.02);
opti.subject_to(tau2(1:end) <= 0.008);

% Lambda parameter - bound on lambda value (different bound for different months)

opti.subject_to(lambda(1:end) >= 0);   % bound on lambda value
opti.subject_to(lambda(1:end) == 0.05);

% Initial Conditions constraints on ODEs

opti.subject_to(S(1,1) == data(1,1));
opti.subject_to(I(1,1) == data(2,1));
opti.subject_to(D(1,1) == data(3,1));
opti.subject_to(T1(1,1) == data(4,1));
opti.subject_to(T2(1,1) == data(5,1));
opti.subject_to(H(1,1) == data(6,1));
opti.subject_to(E(1,1) == data(7,1));

% Initial Conditions constraints on Parameters 

% opti.subject_to(alpha(1,1) == 0.2);
% opti.subject_to(beta(1,1) == 10e-3);
% opti.subject_to(gamma(1,1) == 0.3);


opti.subject_to(S + I + D + T1 + T2 + H + E == 1);

opti.subject_to(S >= 0);
opti.subject_to(I >= 0);
opti.subject_to(D >= 0);
opti.subject_to(T1 >= 0);
opti.subject_to(T2 >= 0);
opti.subject_to(H >= 0);
opti.subject_to(E >= 0);

opti.subject_to(alpha >= 2*beta); % bound on contagion parameter between alpha and beta

% my constraints assumptions 

opti.subject_to(tau1 <= tau2);
opti.subject_to(sigma2 <= sigma1);
opti.subject_to(delta1 > delta2);

%% HARD Constraints on the variation of the coefficents 

% Constraints related to the coefficients - Note that coefficients smoothness will be required only in
% 'constant zones' so in time perioods with constant policies
% For subsequent "dd" value, the constraints will be removed 

% 8 different areas will be covered
policy_idx = [1 40 69 116 141 190 242 294 399];   % Values taken from italian policies applied for the pandemic

for ii=1:length(policy_idx)
    policy_dates(ii) = [new_time(policy_idx(ii))];
end

dd = 1; % days for how much the constraint will be removed

% Note that as assumption I ASSUME that only coefficients \alpha and \beta can be influenced by
% the different policies applied form the government

% constraints NO VARIATION 
for kk = 1:length(policy_idx)-1
    for jj = policy_idx(kk)+dd:policy_idx(kk+1)-1

        opti.subject_to( alpha(1,jj + 1) == alpha(1,jj) )
        opti.subject_to( alpha(1,jj + 1) == alpha(1,jj) )

        opti.subject_to( beta(1,jj + 1) == beta(1,jj) )
        opti.subject_to( beta(1,jj + 1) == beta(1,jj) )

    end
end

% constraint on the first part of the simulation
for jj = 1:policy_idx(2)-1
    opti.subject_to( alpha(1,jj + 1) <= alpha(1,jj)*1.01 )
    opti.subject_to( alpha(1,jj + 1) >= alpha(1,jj)*0.99 )

    opti.subject_to( beta(1,jj + 1) <= beta(1,jj)*1.01 )
    opti.subject_to( beta(1,jj + 1) >= beta(1,jj)*0.99 )
end
%% RK45 simulation of the differential equations

f = @(X, coefs)            [ -X(1)*(coefs(1)*X(2) + coefs(2)*X(3)) ;                                                  % dS/dt = S_dot 
                              ( X(1)*(coefs(1)*X(2) + coefs(2)*X(3)) ) - ( (coefs(3) + coefs(11))*X(2) );             % dI/dt = I_dot
                              (X(2)*coefs(3)) - ( X(3)*(coefs(11) + coefs(4) + coefs(5)) );                           % dD/dt = D_dot
                              ( coefs(4)*X(3) ) - ( (coefs(7) + coefs(9) + coefs(6))*X(4) );                          % dT1/dt = T1_dot
                              ( coefs(5)*X(3) ) - ( (coefs(8) + coefs(10))*X(5) ) + ( coefs(6)*X(4) );                % dT2/dt = T2_dot
                              ( (X(2) + X(3))*coefs(11) ) + ( X(4)*coefs(7) ) + ( X(5)*coefs(8) ) ;                   % dH/dt = H_dot
                              ( X(4)*coefs(9) ) + ( X(5)*coefs(10) )];
                             
dt = 1;

for k = 1:N-1 

% Runge-Kutta 4 integration

   k1 = f(X(:,k),         coefs(:,k));
   k2 = f(X(:,k)+dt/2*k1, coefs(:,k));
   k3 = f(X(:,k)+dt/2*k2, coefs(:,k));
   k4 = f(X(:,k)+dt*k3,   coefs(:,k));
   x_next = X(:,k) + dt/6*(k1+2*k2+2*k3+k4);
   opti.subject_to(X(:,k+1) == x_next);   % close the gaps

end

%% Actual Optimization Simulation

data_obj = horzcat(shortData(:)); 
X_obj = horzcat(X(:));
weightVec = [0.5 0.25 0.75 1 1 0.25 0.1];
Qvec = repmat(weightVec',N,1)';

maxData =( ones(N,1)*max(shortData,[],2)' )';
maxData = horzcat(maxData(:));

coefs_matrix_obj = [    sum( ( (diff(coefs(3,:)))./0.5).^2 );
                        sum( ( (diff(coefs(4,:)))./0.002).^2 );
                        sum( ( (diff(coefs(5,:)))./0.002).^2 );
                        sum( ( (diff(coefs(6,:)))./0.6).^2 );
                        sum( ( (diff(coefs(7,:)))./0.03).^2 );
                        sum( ( (diff(coefs(8,:)))./0.03).^2 );
                        sum( ( (diff(coefs(9,:)))./0.03).^2 );
                        sum( ( (diff(coefs(10,:)))./0.05).^2 );
                        sum( ( (diff(coefs(11,:)))./0.5 ).^2);     
                    ];

cost_matrix_obj =  coefs_matrix_obj(1) + coefs_matrix_obj(2) + coefs_matrix_obj(3) + coefs_matrix_obj(4) + coefs_matrix_obj(5) +...
                   coefs_matrix_obj(6) + coefs_matrix_obj(7) +  coefs_matrix_obj(8) + coefs_matrix_obj(9);

% normalization in cost function - max data
obj = Qvec*((data_obj - X_obj)./maxData).^2 + cost_matrix_obj*1500;
  
opti.minimize(obj);
p_opts = struct('expand', false);
s_opts = struct('max_iter',1e4, 'tol',1e-4,'constr_viol_tol',1e-4,'compl_inf_tol',1e-4,'linear_solver','MA27'); %'MA27' 'MA57' 'MA77' 'MA86' 'MA97'
opti.solver('ipopt',p_opts, s_opts); % set numerical backend
sol = opti.solve();   % actual solver

%% Gathering of all data

% Table of the coefficients

column_names = {'alpha', 'beta', 'gamma', 'delta1', 'delta2', 'epsi', 'sigma1', 'sigma2', 'tau1', 'tau2', 'lambda'};
opti_coefficients = array2table(opti.debug.value(coefs)', 'VariableNames', column_names);

% Trends of every group SIDTTHE
column_names = {'S', 'I', 'D', 'T1', 'T2', 'H', 'E'};
SIDTTHE_trends = array2table(opti.debug.value(X)', 'VariableNames', column_names);

% check if the coefficients reach the bounds

alpha_bound = any(opti_coefficients.alpha == 0.6)
beta_bound = any(opti_coefficients.beta == 0.015)
gamma_bound = any(opti_coefficients.gamma == 0.5)
delta1_bound = any(opti_coefficients.delta1 == 0.002)
delta2_bound = any(opti_coefficients.delta2 == 5e-4)
tau1_bound = any(opti_coefficients.tau1 == 0.03)
tau2_bound = any(opti_coefficients.tau2 == 0.05)
sigma1_bound = any(opti_coefficients.sigma1 == 0.03)
sigma2_bound = any(opti_coefficients.sigma2 == 0.03)
lambda_bound = any(opti_coefficients.lambda == 0.2)

% Saving of the simulation results 

Opti_results = {SIDTTHE_trends, 'Group Trends'; opti_coefficients, 'Coefficients Optimised'; new_time, 'Time for Plotting'};
filename = 'EKF_coefs.mat';
save(filename, 'Opti_results');

%% Coefficients trend plot
close all

% Susceptible - S 
  
figure(1)
plot(new_time, DataArray{1,:}, LineWidth=1.5)
hold on
plot(new_time, SIDTTHE_trends.S, LineWidth=1.5)
ylabel('$\char"0023$ of cases','Interpreter','latex')
title('\textbf{Curve Fitting - Susceptible}','Interpreter','latex')
grid on
legend('Real Data', 'ODE Simulated Data','Interpreter','latex', 'Location','northeast')
xlim([new_time(1), new_time(end)])
set(gca, 'TickLabelInterpreter', 'Latex')

% Infected - I
figure(2)
plot(new_time, DataArray{2,:}, LineWidth=1.5)
hold on
plot(new_time, SIDTTHE_trends.I, LineWidth=1.5)
ylabel('$\char"0023$ of cases','Interpreter','latex')
title('\textbf{Curve Fitting - Infected}','Interpreter','latex')
grid on
legend('Real Data', 'ODE Simulated Data','Interpreter','latex', 'Location','northeast')
xlim([new_time(1), new_time(end)])
set(gca, 'TickLabelInterpreter', 'Latex')


% Diagnosed - D
figure(3)
plot(new_time, DataArray{3,:}, LineWidth=1.5)
hold on
plot(new_time, SIDTTHE_trends.D, LineWidth=1.5)
ylabel('$\char"0023$ of cases','Interpreter','latex')
title('\textbf{Curve Fitting - Diagnosed}','Interpreter','latex')
grid on
legend('Real Data', 'ODE Simulated Data','Interpreter','latex', 'Location','northeast')
xlim([new_time(1), new_time(end)])
set(gca, 'TickLabelInterpreter', 'Latex')


% Lightly Threatned - T1
figure(4)
plot(new_time, DataArray{4,:}, LineWidth=1.5)
hold on
plot(new_time, SIDTTHE_trends.T1, LineWidth=1.5)
ylabel('$\char"0023$ of cases','Interpreter','latex')
title('\textbf{Curve Fitting - Hospitalised }','Interpreter','latex')
grid on
legend('Real Data', 'ODE Simulated Data','Interpreter','latex', 'Location','northeast')
xlim([new_time(1), new_time(end)])
set(gca, 'TickLabelInterpreter', 'Latex')


% ICU / Heavily Threatned - T2
figure(5)
plot(new_time, DataArray{5,:}, LineWidth=1.5)
hold on
plot(new_time, SIDTTHE_trends.T2, LineWidth=1.5)
ylabel('$\char"0023$ of cases','Interpreter','latex')
title('\textbf{Curve Fitting - ICUs / Heavily Threatned }','Interpreter','latex')
grid on
legend('Real Data', 'ODE Simulated Data','Interpreter','latex', 'Location','northeast')
xlim([new_time(1), new_time(end)])
set(gca, 'TickLabelInterpreter', 'Latex')


% Healed - H
figure(6)
plot(new_time, DataArray{6,:}, LineWidth=1.5)
hold on
plot(new_time, SIDTTHE_trends.H, LineWidth=1.5)
ylabel('$\char"0023$ of cases','Interpreter','latex')
title('\textbf{Curve Fitting - Healed}','Interpreter','latex')
grid on
legend('Real Data', 'ODE Simulated Data','Interpreter','latex', 'Location','southeast')
xlim([new_time(1), new_time(end)])
set(gca, 'TickLabelInterpreter', 'Latex')


% Deceased - E
figure(7)
plot(new_time, DataArray{7,:}, LineWidth=1.5)
hold on
plot(new_time, SIDTTHE_trends.E, LineWidth=1.5)
ylabel('$\char"0023$ of cases','Interpreter','latex')
title('\textbf{Curve Fitting - Expired}','Interpreter','latex')
grid on
legend('Real Data', 'ODE Simulated Data','Interpreter','latex', 'Location','southeast')
xlim([new_time(1), new_time(end)])
set(gca, 'TickLabelInterpreter', 'Latex')


% Coefficients alpha
figure(7)
plot(new_time, opti_coefficients.alpha, LineWidth=1.5)
ylabel('Coefficients Values','Interpreter','latex')
title('\textbf{Coefficient $\alpha$}','Interpreter','latex')
grid on
legend('$\alpha$','Interpreter','latex', 'Location','northeast')
xlim([new_time(1), new_time(end)])
set(gca, 'TickLabelInterpreter', 'Latex')


% Coefficients beta
figure(8)
plot(new_time, opti_coefficients.beta, LineWidth=1.5)
ylabel('Coefficients Values','Interpreter','latex')
title('\textbf{Coefficient $\beta$}','Interpreter','latex')
grid on
legend('$\beta$','Interpreter','latex', 'Location','northeast')
xlim([new_time(1), new_time(end)])
set(gca, 'TickLabelInterpreter', 'Latex')


% Coefficient gamma
figure(9)
plot(new_time, opti_coefficients.gamma, LineWidth=1.5)
ylabel('Coefficients Values','Interpreter','latex')
title('\textbf{Coefficient $\gamma$}','Interpreter','latex')
grid on
legend('$\gamma$','Interpreter','latex', 'Location','northeast')
xlim([new_time(1), new_time(end)])
set(gca, 'TickLabelInterpreter', 'Latex')


% Coefficients delta1 and delta2
figure(10)
plot(new_time, opti_coefficients.delta1, LineWidth=1.5)
hold on
plot(new_time, opti_coefficients.delta2, LineWidth=1.5)
ylabel('Coefficients Values','Interpreter','latex')
title('\textbf{Coefficients $\delta_1$ and $\delta_2$}','Interpreter','latex')
grid on
legend('$\delta_1$', '$\delta_2$','Interpreter','latex', 'Location','northeast')
xlim([new_time(1), new_time(end)])
set(gca, 'TickLabelInterpreter', 'Latex')


% Coefficients sigma1 and sigma2
figure(12)
plot(new_time, opti_coefficients.sigma1, LineWidth=1.5)
hold on
plot(new_time, opti_coefficients.sigma2, LineWidth=1.5)
ylabel('Coefficients Values','Interpreter','latex')
title('\textbf{Coefficients $\sigma_1$ and $\sigma_2$}','Interpreter','latex')
grid on
legend('$\sigma_1$','$\sigma_2$','Interpreter','latex', 'Location','northeast')
xlim([new_time(1), new_time(end)])
set(gca, 'TickLabelInterpreter', 'Latex')


% Coefficients tau1 and tau2
 figure(13)
plot(new_time, opti_coefficients.tau1, LineWidth=1.5)
hold on
plot(new_time, opti_coefficients.tau2, LineWidth=1.5)
ylabel('Coefficients Values','Interpreter','latex')
title('\textbf{Coefficients $\tau_1$ and $\tau_2$}','Interpreter','latex')
grid on
legend('$\tau_1$','$\tau_2$','Interpreter','latex', 'Location','northeast')
xlim([new_time(1), new_time(end)])
set(gca, 'TickLabelInterpreter', 'Latex')


% Coefficients lambda
figure(14)
plot(new_time, opti_coefficients.lambda, LineWidth=1.5)
ylabel('Coefficients Values','Interpreter','latex')
title('\textbf{Coefficient $\lambda$}','Interpreter','latex')
grid on
legend('$\lambda$', 'Interpreter','latex', 'Location','northeast')
xlim([new_time(1), new_time(end)])
set(gca, 'TickLabelInterpreter', 'Latex')
