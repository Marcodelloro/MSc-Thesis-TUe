clear all
close all
clc
addpath('/Users/marcodelloro/Downloads/casadi-3.6.3-osx64-matlab2018b')
import casadi.*;
opti = casadi.Opti();

%------------------- Fit the SIDTTHE model with the data - CasADi Code --------------------%


%% VERSION 1 - Epsilon presence and exchanges between T1 and T2

load('/Users/marcodelloro/Desktop/Thesis/MSc-Thesis-TUe/Data_Colection/SIDTTHE_data.mat');

% Normalisation of loaded data so they add up to 1

date = SIDTTHE_data{1,1}.date;
Npop = 59240329; % Total Population of Italy

I_data = SIDTTHE_data{1,1}.data / Npop;     % I group - Infected Population **
D_data = SIDTTHE_data{2,1}.data / Npop;     % D group - Diagnosed Population
T1_data = SIDTTHE_data{3,1}.data / Npop;    % T1 group - Hospitalised (Lightly Threatned) Population
T2_data = SIDTTHE_data{4,1}.data / Npop;    % T2 group - ICUs (Severly Threatned) Population
H_data = SIDTTHE_data{6,1}.data / Npop;     % T2 group - Healed Population **
E_data = SIDTTHE_data{5,1}.data / Npop;     % E group - Expired (Deceased) Population

S_data = ones(length(I_data),1) - (I_data + D_data + T1_data + T2_data + H_data + E_data); 

% ** the groups I and H are based on strong assumptions and dark numbers estimates

DataArray = {S_data, I_data, D_data, T1_data, T2_data, H_data, E_data};

% Upsampling of the dataset using Interp1, so that we can handle better the divided day by day 

% New time vector (daily data)
new_time = linspace(SIDTTHE_data{1, 1}.date(1), SIDTTHE_data{1, 1}.date(end), 399); % 399 data points (daily data)

% Iterate over each cell in DataArray
for ii = 1:numel(DataArray)

    original_data = DataArray{1, ii};

    upsampled_DataArray{ii,1} = interp1(SIDTTHE_data{1, 1}.date, original_data, new_time);

end

%% Setting of the optimization variables

N = 399; % if we are using the DAILY time discretization, then N = 399
         % else is N = 57 (weekly)

alpha = opti.variable(1,N);   % alpha - coefficient to go from S to I 
beta = opti.variable(1,N);    % beta - coefficient to go from S to D
gamma = opti.variable(1,N);   % gamma - coefficient to go from I to D
delta1 = opti.variable(1,N);  % delta 1 - coefficient to go from D to T1
delta2 = opti.variable(1,N);  % delta 2 - coefficient to go from D to T2 
epsi1 = opti.variable(1,N);   % epsilon 1 - coefficient to go from T1 to T2
epsi2 = opti.variable(1,N);   % epsilon 2 - coefficient to go from T2 to T1
sigma1 = opti.variable(1,N);  % sigma 1 - coefficient to go from T1 to H
sigma2 = opti.variable(1,N);  % sigma 2 - coefficient to go from T2 to H
tau1 = opti.variable(1,N);    % tau 1 - coefficient to go from T1 to E
tau2 = opti.variable(1,N);    % tau 2 - coefficient to go from T2 to E
lambda = opti.variable(1,N);  % lambda - coefficient to go from I and D to H

coefs = [alpha; beta; gamma; delta1; delta2; epsi1; epsi2; sigma1; sigma2; tau1; tau2; lambda];

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
time_intervals = 399;

opti.set_initial(S,ones(time_intervals,1) * S_data(1));
opti.set_initial(I,ones(time_intervals,1) * I_data(1));
opti.set_initial(D,ones(time_intervals,1) * D_data(1));
opti.set_initial(T1,ones(time_intervals,1)* T1_data(1));
opti.set_initial(T2,ones(time_intervals,1)* T2_data(1));
opti.set_initial(H,ones(time_intervals,1) * H_data(1));
opti.set_initial(E,ones(time_intervals,1) * E_data(1));

opti.set_initial(alpha,ones(time_intervals,1) * 0.513);          % All coeff. are obtained from the mean of the values provided by G.Giordano in original SIDARTHE paper
opti.set_initial(beta,ones(time_intervals,1) * 0.011); 
opti.set_initial(gamma,ones(time_intervals,1) * 0.2405);
opti.set_initial(delta1,ones(time_intervals,1) * 0.01);
opti.set_initial(delta2,ones(time_intervals,1) * 0.008);
% opti.set_initial(epsi1,H.data(1));                             % No initial values for this two, will be totally estimated 
% opti.set_initial(epsi2,E.data(1)); 
opti.set_initial(sigma1,ones(time_intervals,1) * 0.03);
opti.set_initial(sigma2,ones(time_intervals,1) * 0.065);
opti.set_initial(tau1,ones(time_intervals,1) * 0.015);
opti.set_initial(tau2,ones(time_intervals,1) * 0.005);
opti.set_initial(lambda,ones(time_intervals,1) * 0.0596);

% Constraints for the parameters

opti.subject_to(alpha >= 0);
opti.subject_to(alpha <= 0.6);   % bound on alfa value

opti.subject_to(beta >= 0);
opti.subject_to(beta <= 0.015); % bound on beta value

opti.subject_to(gamma >= 0);    % bound on gamma value
opti.subject_to(gamma <= 0.5);

opti.subject_to(delta1 >= 0.001);   % bound on delta1 value
opti.subject_to(delta1 <= 0.02);

opti.subject_to(delta2 >= 0);   % bound on delta2 value
opti.subject_to(delta2 <= 0.03);

opti.subject_to(epsi1 >= 0);    % bound on epsi1 value
opti.subject_to(epsi1 <= 0.5);

opti.subject_to(epsi2 >= 0);    % bound on epsi2 value
opti.subject_to(epsi2 <= 0.6);

opti.subject_to(sigma1 >= 0);   % bound on sigma1 value
opti.subject_to(sigma1 <= 0.03);

opti.subject_to(sigma2 >= 0);   % bound on sigma2 value
opti.subject_to(sigma2 <= 0.03);

opti.subject_to(tau1 >= 0);     % bound on tau1 value
opti.subject_to(tau1 <= 0.03);

opti.subject_to(tau2 >= 0);     % bound on tau2 value
% opti.subject_to(tau2 <= 0.02);
opti.subject_to(tau2 <= 0.05);

opti.subject_to(lambda >= 0);   % bound on lambda value
opti.subject_to(lambda <= 0.2);

opti.subject_to(S + I + D + T1 + T2 + H + E == 1);

opti.subject_to(S >= 0);
opti.subject_to(I >= 0);
opti.subject_to(D >= 0);
opti.subject_to(T1 >= 0);
opti.subject_to(T2 >= 0);
opti.subject_to(H >= 0);
opti.subject_to(E >= 0);

opti.subject_to(alpha >= 3*beta); % bound on contagion parameter between alpha and beta

% my constraints assumptions 

opti.subject_to(epsi2 <= epsi1);
opti.subject_to(tau1 >= tau2);
opti.subject_to(sigma2 >= sigma1);
opti.subject_to(delta1 > 1.5 * delta2);

% constraints on the variation of the coefficents

for ii = 1:N-1

    opti.subject_to( alpha(1,ii + 1) <= alpha(1,ii)*1.01 )
    opti.subject_to( alpha(1,ii + 1) >= alpha(1,ii)*0.99 )

    opti.subject_to( beta(1,ii + 1) <= beta(1,ii)*1.05 )
    opti.subject_to( beta(1,ii + 1) >= beta(1,ii)*0.95 )

    opti.subject_to( gamma(1,ii + 1) <= gamma(1,ii)*1.01 )
    opti.subject_to( gamma(1,ii + 1) >= gamma(1,ii)*0.99 )

    opti.subject_to( delta1(1,ii + 1) <= delta1(1,ii)*1.05 )
    opti.subject_to( delta1(1,ii + 1) >= delta1(1,ii)*0.95 )

    opti.subject_to( delta2(1,ii + 1) <= delta2(1,ii)*1.05 )
    opti.subject_to( delta2(1,ii + 1) >= delta2(1,ii)*0.95 )

    opti.subject_to( epsi1(1,ii + 1) <= epsi1(1,ii)*1.05 )
    opti.subject_to( epsi1(1,ii + 1) >= epsi1(1,ii)*0.95 )

    opti.subject_to( epsi2(1,ii + 1) <= epsi2(1,ii)*1.05 )
    opti.subject_to( epsi2(1,ii + 1) >= epsi2(1,ii)*0.95)

    opti.subject_to( sigma1(1,ii + 1) <= sigma1(1,ii)*1.05 )
    opti.subject_to( sigma1(1,ii + 1) >= sigma1(1,ii)*0.95 )

    opti.subject_to( sigma2(1,ii + 1) <= sigma2(1,ii)*1.05 )
    opti.subject_to( sigma2(1,ii + 1) >= sigma2(1,ii)*0.95 )

    opti.subject_to( tau1(1,ii + 1) <= tau1(1,ii)*1.05 )
    opti.subject_to( tau1(1,ii + 1) >= tau1(1,ii)*0.95 )

    opti.subject_to( tau2(1,ii + 1) <= tau2(1,ii)*1.05 )
    opti.subject_to( tau2(1,ii + 1) >= tau2(1,ii)*0.95)

    opti.subject_to( lambda(1,ii + 1) <= lambda(1,ii)*1.05 )
    opti.subject_to( lambda(1,ii + 1) >= lambda(1,ii)*0.95 )

end

%% 
% Specify system dynamics - creation of the handle functions

% Version 1
f = @(X, coefs)   [           -X(1).*(coefs(1)*X(2) + coefs(2)*X(3));                                                 % dS/dt = S_dot 
                              ( X(1)*(coefs(1)*X(2) + coefs(2)*X(3)) ) - ( (coefs(3) + coefs(12))*X(2) );             % dI/dt = I_dot
                              (X(2)*coefs(3)) - ( X(3)*(coefs(12) + coefs(4) + coefs(5)) );                           % dD/dt = D_dot
                              (coefs(4)*X(3)) - ( (coefs(8) + coefs(10) + coefs(6))*X(4) ) + (coefs(7)*X(5)) ;        % dT1/dt = T1_dot
                              (coefs(5)*X(3)) - ( (coefs(9) + coefs(11) + coefs(7))*X(5) ) + (coefs(6)*X(4)) ;        % dT2/dt = T2_dot
                              ( (X(2) + X(3))*coefs(12) ) + (X(4)*coefs(8)) + (X(5)*coefs(9)) ;                       % dH/dt = H_dot
                              (X(4)*coefs(10)) + (X(5)*coefs(11)) ];                                                  % dE/dt = E_dot


dt = 1;
for k=1:N-1      % loop over control intervals
                
% Runge-Kutta 4 integration

   k1 = f(X(:,k),         coefs(:,k));
   k2 = f(X(:,k)+dt/2*k1, coefs(:,k));
   k3 = f(X(:,k)+dt/2*k2, coefs(:,k));
   k4 = f(X(:,k)+dt*k3,   coefs(:,k));
   x_next = X(:,k) + dt/6*(k1+2*k2+2*k3+k4);
   opti.subject_to(X(:,k+1) == x_next); % close the gaps
end

%% Actual Optimization Simulation

data = vertcat(upsampled_DataArray{:});
data_obj = horzcat(data(:));
X_obj = horzcat(X(:));

obj = sum((data_obj - X_obj).^2)*1e3 ;
opti.minimize(obj);

%options = struct('max_iter',1e4,'tol',1e-5,'linear_solver','MA27');
opti.solver('ipopt'); % set numerical backend
sol = opti.solve();   % actual solver

% Gathering of all data

% Table of the coefficients
column_names = {'alpha', 'beta', 'gamma', 'delta1', 'delta2', 'epsi1', 'epsi2', 'sigma1', 'sigma2', 'tau1', 'tau2', 'lambda'};
opti_coefficients = array2table(opti.debug.value(coefs)', 'VariableNames', column_names);

% Trends of every group SIDTTHE
column_names = {'S', 'I', 'D', 'T1', 'T2', 'H', 'E'};
SIDTTHE_trends = array2table(opti.debug.value(X)', 'VariableNames', column_names);

% check if the coefficients reach the bounds

alpha_bound = any(opti_coefficients.alpha == 0.6)
beta_bound = any(opti_coefficients.beta == 0.015)
gamma_bound = any(opti_coefficients.gamma == 0.5)
delta1_bound = any(opti_coefficients.delta1 == 0.02)
delta2_bound = any(opti_coefficients.delta2 == 0.03)
tau1_bound = any(opti_coefficients.tau1 == 0.03)
tau2_bound = any(opti_coefficients.tau2 == 0.05)
sigma1_bound = any(opti_coefficients.sigma1 == 0.03)
sigma2_bound = any(opti_coefficients.sigma2 == 0.03)
lambda_bound = any(opti_coefficients.lambda == 0.2)

% Saving of the simulation results 

Opti_results = {SIDTTHE_trends, 'Group Trends'; opti_coefficients, 'Coefficients Optimised'; new_time, 'Time for Plotting'};
filename = 'Optimization_Results.mat';
save(filename, 'Opti_results');

%% Comparison Plot

% Susceptible - S 
figure(1)
plot(new_time, upsampled_DataArray{1,:}, LineWidth=1.5)
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
plot(new_time, upsampled_DataArray{2,:}, LineWidth=1.5)
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
plot(new_time, upsampled_DataArray{3,:}, LineWidth=1.5)
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
plot(new_time, upsampled_DataArray{4,:}, LineWidth=1.5)
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
plot(new_time, upsampled_DataArray{5,:}, LineWidth=1.5)
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
plot(new_time, upsampled_DataArray{6,:}, LineWidth=1.5)
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
plot(new_time, upsampled_DataArray{7,:}, LineWidth=1.5)
hold on
plot(new_time, SIDTTHE_trends.E, LineWidth=1.5)
ylabel('$\char"0023$ of cases','Interpreter','latex')
title('\textbf{Curve Fitting - Expired}','Interpreter','latex')
grid on
legend('Real Data', 'ODE Simulated Data','Interpreter','latex', 'Location','southeast')
xlim([new_time(1), new_time(end)])
set(gca, 'TickLabelInterpreter', 'Latex')

%% Coefficients trend plot


% Coefficients alpha and beta
figure(8)
plot(new_time, opti_coefficients.alpha, LineWidth=1.5)
hold on
plot(new_time, opti_coefficients.beta, LineWidth=1.5)
ylabel('Coefficients Values','Interpreter','latex')
title('\textbf{Coefficients $\alpha$ and $\beta$}','Interpreter','latex')
grid on
legend('$\alpha$', '$\beta$','Interpreter','latex', 'Location','northeast')
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


% Coefficients epsi1 and epsi2
figure(11)
plot(new_time, opti_coefficients.epsi1, LineWidth=1.5)
hold on
plot(new_time, opti_coefficients.epsi2, LineWidth=1.5)
ylabel('Coefficients Values','Interpreter','latex')
title('\textbf{Coefficients $\epsilon_1$ and $\epsilon_2$}','Interpreter','latex')
grid on
legend('$\epsilon_1$','$\epsilon_2$','Interpreter','latex', 'Location','northeast')
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
figure(13)
plot(new_time, opti_coefficients.lambda, LineWidth=1.5)
ylabel('Coefficients Values','Interpreter','latex')
title('\textbf{Coefficients $\tau_1$ and $\tau_2$}','Interpreter','latex')
grid on
legend('$\lambda$', 'Interpreter','latex', 'Location','northeast')
xlim([new_time(1), new_time(end)])
set(gca, 'TickLabelInterpreter', 'Latex')