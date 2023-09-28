clear all
close all
clc
addpath('/Users/marcodelloro/Downloads/casadi-3.6.3-osx64-matlab2018b')
import casadi.*;
opti = casadi.Opti();

%% Fit the SIDTTHE model with the data - CasADi Code

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

% S_dot = opti.variable(1,N);
% I_dot = opti.variable(1,N);
% D_dot = opti.variable(1,N);
% T1_dot = opti.variable(1,N);
% T2_dot = opti.variable(1,N);
% H_dot = opti.variable(1,N);
% E_dot = opti.variable(1,N);

% Initial conditions for the various parameters
time_intervals = 399;

opti.set_initial(S,ones(time_intervals,1) * S_data(1));
opti.set_initial(I,ones(time_intervals,1) * I_data(1));
opti.set_initial(D,ones(time_intervals,1) * D_data(1));
opti.set_initial(T1,ones(time_intervals,1)* T1_data(1));
opti.set_initial(T2,ones(time_intervals,1)* T2_data(1));
opti.set_initial(H,ones(time_intervals,1) * H_data(1));
opti.set_initial(E,ones(time_intervals,1) * E_data(1));

opti.set_initial(alpha,0.513);          % All coeff. are obtained from the mean of the values provided by G.Giordano in original SIDARTHE paper
opti.set_initial(beta,0.011); 
opti.set_initial(gamma,0.2405);
opti.set_initial(delta1,0.01);
opti.set_initial(delta2,0.008);
% opti.set_initial(epsi1,H.data(1));    % No initial values for this two, will be totally estimated 
% opti.set_initial(epsi2,E.data(1)); 
opti.set_initial(sigma1,0.03);
opti.set_initial(sigma2,0.065);
opti.set_initial(tau1,0.015);
opti.set_initial(tau2,0.005);
opti.set_initial(lambda,0.0596);

% Constraints for the parameters

opti.subject_to(alpha >= 0);
opti.subject_to(alpha <= 0.6);   % bound on alfa value

opti.subject_to(beta >= 0);
opti.subject_to(beta <= 0.015); % bound on beta value

opti.subject_to(gamma >= 0);    % bound on gamma value
opti.subject_to(gamma <= 0.5);

opti.subject_to(delta1 >= 0);   % bound on delta1 value
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
opti.subject_to(tau2 <= 0.02);

opti.subject_to(lambda >= 0);   % bound on lambda value
opti.subject_to(lambda <= 0.2);

opti.subject_to(S + I + D + T1 + T2 + H + E == 1);

opti.subject_to(alpha >= 3*beta); % bound on contagion parameter between alpha and beta

% my constraints assumptions 

opti.subject_to(epsi2 >= epsi1);
opti.subject_to(tau1 >= tau2);
opti.subject_to(sigma2 >= sigma1);

%% 
% Specify system dynamics - creation of the handle functions

f = @(X, coefs)   [           -X(1).*(coefs(1)*X(2) + coefs(2)*X(3));                                                 % dS/dt = S_dot 
                              ( X(1)*(coefs(1)*X(2) + coefs(2)*X(3)) ) - ( (coefs(3) + coefs(12))*X(2) );             % dI/dt = I_dot
                              (X(2)*coefs(3)) - ( X(3)*(coefs(12) + coefs(4) + coefs(5)) );                           % dD/dt = D_dot
                              (coefs(4)*X(3)) - ( (coefs(8) + coefs(10) + coefs(6))*X(4) ) + (coefs(7)*X(5)) ;         % dT1/dt = T1_dot
                              (coefs(5)*X(3)) - ( (coefs(9) + coefs(11) + coefs(7))*X(5) ) + (coefs(6)*X(4)) ;        % dT2/dt = T2_dot
                              ( (X(2) + X(3))*coefs(12) ) + (X(4)*coefs(8)) + (X(5)*coefs(9)) ;                       % dH/dt = H_dot
                              (X(4)*coefs(10)) + (X(5)*coefs(11)) ];                                                  % dE/dt = E_dot



% % Specify system dynamics - creation of the handle functions
% 
% f = @(X, coefs)                                                     [-S.*(alpha*I + beta*D)];                                         % dS/dt = S_dot
% 
% I_dot = @(S, I, D, alfa, beta, gamma, lambda)                    [( S*(alpha*I + beta*D) ) - ( (gamma + lambda)*I )];             % dI/dt = I_dot
% 
% D_dot = @(I, D, gamma, lambda, delta1, delta2)                   [(I*gamma) - ( D*(lambda + delta1 + delta2) )];                  % dD/dt = D_dot
% 
% T1_dot = @(D, T1, T2, delta1, sigma1, epsi1, tau1, epsi2)        [(delta1*D) - ( (sigma1 + tau1 + epsi1)*T1 ) + (epsi2*T2) ];     % dT1/dt = T1_dot
% 
% T2_dot = @(D, T1, T2, delta2, sigma2, epsi1, tau2, epsi2)        [(delta2*D) - ( (sigma2 + tau2 + epsi2)*T2 ) + (epsi1*T1) ];     % dT2/dt = T2_dot
% 
% H_dot = @(I, D, T1, T2, lambda, sigma1, sigma2)                  [ ( (I + D)*lambda ) + (T1*sigma1) + (T2*sigma2) ];              % dH/dt = H_dot
% 
% E_dot = @(T1, T2, tau1, tau2 )                                   [ (T1*tau1) + (T2*tau2) ];                                       % dE/dt = E_dot


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
obj = sum((data_obj - X_obj).^2);
opti.minimize(obj);

opti.solver('ipopt'); % set numerical backend
sol = opti.solve();   % actual solver

%% Gathering of all data

% Table of the coefficients
column_names = {'alpha', 'beta', 'gamma', 'delta1', 'delta2', 'epsi1', 'epsi2', 'sigma1', 'sigma2', 'tau1', 'tau2', 'lambda'};
opti_coefficients = array2table(opti.debug.value(coefs)', 'VariableNames', column_names);

% Trends of every group SIDTTHE
column_names = {'S', 'I', 'D', 'T1', 'T2', 'H', 'E'};
SIDTTHE_trends = array2table(opti.debug.value(X)', 'VariableNames', column_names);


%% Comparison Plot

% Susceptible - S 
figure()
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
figure()
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
figure()
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
figure()
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
figure()
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
figure()
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
figure()
plot(new_time, upsampled_DataArray{7,:}, LineWidth=1.5)
hold on
plot(new_time, SIDTTHE_trends.E, LineWidth=1.5)
ylabel('$\char"0023$ of cases','Interpreter','latex')
title('\textbf{Curve Fitting - Expired}','Interpreter','latex')
grid on
legend('Real Data', 'ODE Simulated Data','Interpreter','latex', 'Location','southeast')
xlim([new_time(1), new_time(end)])
set(gca, 'TickLabelInterpreter', 'Latex')