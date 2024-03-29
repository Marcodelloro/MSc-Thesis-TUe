 %% Version with VACCINES - Closer to the final model

% This version of the model should contain also the data coming from
% vaccinated individuals , the model is SIDDTHE-V model - MOD with
% POPULATION COMING FROM H and S 

clc
clear all

load('/Users/marcodelloro/Desktop/Thesis/MSc-Thesis-TUe/Data_Collection/SIDTTHE_data.mat');
load('/Users/marcodelloro/Desktop/Thesis/MSc-Thesis-TUe/Data_Collection/Tests_data.mat');
load('/Users/marcodelloro/Desktop/Thesis/MSc-Thesis-TUe/Data_Collection/vaxData.mat'); 

addpath('/Users/marcodelloro/Downloads/casadi-3.6.3-osx64-matlab2018b')
import casadi.*;
opti = casadi.Opti();

% Normalisation of loaded data so they add up to 1

date = SIDTTHE_data{1,1}.date;
Npop = 59240329; % Total Population of Italy

I_data = SIDTTHE_data{1,1}.data / Npop;     % I group - Infected Population **
D_data = SIDTTHE_data{2,1}.data / Npop;     % D group - Diagnosed Population
T1_data = SIDTTHE_data{3,1}.data / Npop;    % T1 group - Hospitalised (Lightly Threatned) Population
T2_data = SIDTTHE_data{4,1}.data / Npop;    % T2 group - ICUs (Severly Threatned) Population
H_data = (SIDTTHE_data{6,1}.data - vaxData.sum_dpi') / Npop;     % T2 group - Healed Population **
E_data = SIDTTHE_data{5,1}.data  / Npop;     % E group - Expired (Deceased) Population
V_data = (vaxData.sum_d2 + vaxData.sum_dpi)' / Npop;
S_data = ones(length(I_data),1)' - (I_data + D_data + T1_data + T2_data + H_data + E_data + V_data ); 

% ** the groups I and H are based on strong assumptions and dark numbers estimates

DataArray = {S_data ; I_data ; D_data ; T1_data ; T2_data ; H_data ; E_data ; V_data };
N = 399;

data = cell2mat(DataArray); % to use non-scored data

new_time = date;

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
phi = opti.variable(1,N);  % phi - coefficient to go from S to V
kappa = opti.variable(1,N); % phi - coefficient to go from H to V

coefs = [alpha; beta; gamma; delta1; delta2; epsi; sigma1; sigma2; tau1; tau2; lambda; phi; kappa];

% Lifting Variables

S = opti.variable(1,N); 
I = opti.variable(1,N);
D = opti.variable(1,N);
T1 = opti.variable(1,N);
T2 = opti.variable(1,N);
H = opti.variable(1,N);
E = opti.variable(1,N);
V = opti.variable(1,N);

X = [S; I; D; T1; T2; H; E; V]; % state trajectory

% Initial conditions for the various parameters

opti.set_initial(S,ones(N,1) * data(1,1));
opti.set_initial(I,ones(N,1) * data(2,1));
opti.set_initial(D,ones(N,1) * data(3,1));
opti.set_initial(T1,ones(N,1)* data(4,1));
opti.set_initial(T2,ones(N,1)* data(5,1));
opti.set_initial(H,ones(N,1) * data(6,1));
opti.set_initial(E,ones(N,1) * data(7,1));
opti.set_initial(V,ones(N,1) * data(8,1));

opti.set_initial(alpha,ones(N,1) * 0.513);          % All coeff. are obtained from the mean of the values provided by G.Giordano in original SIDARTHE paper
% opti.set_initial(beta,ones(N,1) * 0.011);
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

% opti.subject_to(beta(1:end) >= 0);
% opti.subject_to(beta(1:end) <= 0.015);

opti.subject_to(beta(1:end) == 1e-6);

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
opti.subject_to(lambda(1:end) <= 0.2);

% phi parameter - bound on phi value (different bound for different months)

opti.subject_to(phi(1:end) >= 0);   % bound on phi value
opti.subject_to(kappa(1:end) >= 0);   % bound on phi value


% Initial Conditions constraints on ODEs

opti.subject_to(S(1,1) == data(1,1));
opti.subject_to(I(1,1) == data(2,1));
opti.subject_to(D(1,1) == data(3,1));
opti.subject_to(T1(1,1) == data(4,1));
opti.subject_to(T2(1,1) == data(5,1));
opti.subject_to(H(1,1) == data(6,1));
opti.subject_to(E(1,1) == data(7,1));
opti.subject_to(V(1,1) == data(8,1));

% Initial Conditions constraints on Parameters 

% opti.subject_to(alpha(1,1) == 0.2);
% opti.subject_to(beta(1,1) == 10e-3);
% opti.subject_to(gamma(1,1) == 0.3);


opti.subject_to(S + I + D + T1 + T2 + H + E + V == 1);

opti.subject_to(S >= 0);
opti.subject_to(I >= 0);
opti.subject_to(D >= 0);
opti.subject_to(T1 >= 0);
opti.subject_to(T2 >= 0);
opti.subject_to(H >= 0);
opti.subject_to(E >= 0);
opti.subject_to(V >= 0);

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

%---------------------------------------------------------------------
% 
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
for jj = 1:policy_idx(2)-1
    opti.subject_to( alpha(1,jj + 1) <= alpha(1,jj)*1.01 )
    opti.subject_to( alpha(1,jj + 1) >= alpha(1,jj)*0.99 )

    % opti.subject_to( beta(1,jj + 1) <= beta(1,jj)*1.01 )
    % opti.subject_to( beta(1,jj + 1) >= beta(1,jj)*0.99 )
end

 % opti.subject_to( alpha(1,policy_idx(2)+dd) <= alpha(1,policy_idx(2))*0.9 )
 % opti.subject_to( beta(1,policy_idx(2)+dd) <= beta(1,policy_idx(2))*0.9 )
 % opti.subject_to( beta(1,policy_idx(3)+dd) <= beta(1,policy_idx(3)) )
 % opti.subject_to( alpha(1,policy_idx(7)+dd) >= alpha(1,policy_idx(7))*1.1)
 % opti.subject_to( beta(1,policy_idx(7)+dd) > beta(1,policy_idx(7))*1.1 )
 % opti.subject_to( beta(1,policy_idx(6)+dd) <= beta(1,policy_idx(6))*0.9 )
 % opti.subject_to( alpha(1,policy_idx(6)+dd) <= alpha(1,policy_idx(6))*0.9 )
 % opti.subject_to( beta(1,policy_idx(8)+dd) >= beta(1,policy_idx(8))*1.1 )
 % opti.subject_to( beta(1,policy_idx(4)+dd) <= beta(1,policy_idx(4))*0.9 )
 % opti.subject_to( alpha(1,policy_idx(4)+dd) <= alpha(1,policy_idx(4))*0.9 )

%% RK45 simulation of the differential equations

% Specify system dynamics - creation of the handle functions

f = @(X, coefs, q)            [ -X(1)*(coefs(1)*X(2) + coefs(2)*X(3)) - X(1)*coefs(12) ;                                                 % dS/dt = S_dot 
                              ( X(1)*(coefs(1)*X(2) + coefs(2)*X(3)) ) - ( (coefs(3) + coefs(11))*X(2) );             % dI/dt = I_dot
                              (X(2)*coefs(3)) - ( X(3)*(coefs(11) + coefs(4) + coefs(5)) );                           % dD/dt = D_dot
                              ( coefs(4)*X(3) ) - ( (coefs(7) + coefs(9) + coefs(6))*X(4) );                            % dT1/dt = T1_dot
                              ( coefs(5)*X(3) ) - ( (coefs(8) + coefs(10))*X(5) ) + ( coefs(6)*X(4) );                  % dT2/dt = T2_dot
                              ( (X(2) + X(3))*coefs(11) ) + ( X(4)*coefs(7) ) + ( X(5)*coefs(8) ) - X(6)*coefs(13) ;                       % dH/dt = H_dot
                              ( X(4)*coefs(9) ) + ( X(5)*coefs(10) );
                                X(1)*coefs(12) + X(6)*coefs(13) ];
                             

dt = 1;

for k = 1:N-1 

% Runge-Kutta 4 integration

   k1 = f(X(:,k),         coefs(:,k));
   k2 = f(X(:,k)+dt/2*k1, coefs(:,k));
   k3 = f(X(:,k)+dt/2*k2, coefs(:,k));
   k4 = f(X(:,k)+dt*k3,   coefs(:,k));
   x_next = X(:,k) + dt/6*(k1+2*k2+2*k3+k4) + [exp(-1/180)*X(8); zeros(7,1)];
   opti.subject_to(X(:,k+1) == x_next); % close the gaps

end


%% Actual Optimization Simulation

data_obj = horzcat(data(:)); 
X_obj = horzcat(X(:));
weightVec = [0.5 0.25 0.75 1 1 0.25 0.1 0.75];
Qvec = repmat(weightVec',399,1)';

maxData =( ones(399,1)*max(data,[],2)' )';
maxData = horzcat(maxData(:));

meanData =( ones(399,1)*mean(data,2)' )';
meanData = horzcat(meanData(:));


% building the actual COST FUNCTIONS to ensure smoothing of unconstrained coefs

coefs_matrix_obj = [    sum( ( (diff(coefs(3,:)))./0.5).^2 );
                        sum( ( (diff(coefs(4,:)))./0.002).^2 );
                        sum( ( (diff(coefs(5,:)))./0.002).^2 );
                        sum( ( (diff(coefs(6,:)))./0.6).^2 );
                        sum( ( (diff(coefs(7,:)))./0.03).^2 );
                        sum( ( (diff(coefs(8,:)))./0.03).^2 );
                        sum( ( (diff(coefs(9,:)))./0.03).^2 );
                        sum( ( (diff(coefs(10,:)))./0.05).^2 );
                        sum( ( (diff(coefs(11,:)))./0.5 ).^2);     
                        sum( ( (diff(coefs(12,:)))./0.5 ).^2);  
                        sum( ( (diff(coefs(13,:)))./0.5 ).^2)     ];

cost_matrix_obj =  coefs_matrix_obj(1) + coefs_matrix_obj(2) + coefs_matrix_obj(3) + coefs_matrix_obj(4) + coefs_matrix_obj(5) +...
                   coefs_matrix_obj(6) + coefs_matrix_obj(7) +  coefs_matrix_obj(8) + coefs_matrix_obj(9) +  coefs_matrix_obj(10) + coefs_matrix_obj(11);


integral_cost1= 0; % initialization of the integral cost for Hosipitalised
integral_cost2= 0; % initialization of the integral cost for ICUs

for ii = 1:N
    step_cost1 = (data(4,ii) - X(4,ii)).^2;
    step_cost2 = (data(5,ii) - X(5,ii)).^2;
    integral_cost1 = integral_cost1 + step_cost1;
    integral_cost2 = integral_cost2 + step_cost2;
end

% Definitions of some weights for the optimization

w1 = 6;     % weight for the difference between data and model
w2 = 10;     % weight for the integral cost 1 
w3 = 1;     % weight for the integral cost 2
w4 = 20;   % weight on coefficient smoothing 

% normalization in cost function - mean data
%  obj = sum(((data_obj - X_obj)./meanData).^2)*0.01 + cost_matrix_obj*100;

% normalization in cost function - max data

obj = Qvec*((data_obj - X_obj)./maxData).^2 + cost_matrix_obj*100;

opti.minimize(obj);
p_opts = struct('expand', false);
s_opts = struct('max_iter',1e4, 'tol',1e-4,'constr_viol_tol',1e-4,'compl_inf_tol',1e-4,'linear_solver','MA97'); %'MA27' 'MA57' 'MA77' 'MA86' 'MA97'
opti.solver('ipopt',p_opts, s_opts); % set numerical backend
sol = opti.solve();   % actual solver

%% Gathering of all data

% Table of the coefficients

column_names = {'alpha', 'beta', 'gamma', 'delta1', 'delta2', 'epsi', 'sigma1', 'sigma2', 'tau1', 'tau2', 'lambda','phi','kappa'};
opti_coefficients = array2table(opti.debug.value(coefs)', 'VariableNames', column_names);

% Trends of every group SIDTTHE
column_names = {'S', 'I', 'D', 'T1', 'T2', 'H', 'E', 'V'};
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
filename = 'Optimization_Results.mat';
save(filename, 'Opti_results');


%% Comparison Plot

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

% Vaccinated - V
figure(8)
plot(new_time, DataArray{8,:}, LineWidth=1.5)
hold on
plot(new_time, SIDTTHE_trends.V, LineWidth=1.5)
ylabel('$\char"0023$ of cases','Interpreter','latex')
title('\textbf{Curve Fitting - Vaccinated}','Interpreter','latex')
grid on
legend('Real Data', 'ODE Simulated Data','Interpreter','latex', 'Location','southeast')
xlim([new_time(1), new_time(end)])
set(gca, 'TickLabelInterpreter', 'Latex')

%% Coefficients trend plot
close all

% Coefficients alpha
figure(7)
plot(new_time, opti_coefficients.alpha, LineWidth=1.5)
ylabel('Coefficients Values','Interpreter','latex')
title('\textbf{Coefficient $\alpha$}','Interpreter','latex')
grid on
legend('$\beta$','Interpreter','latex', 'Location','northeast')
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
ylim([0, max(opti_coefficients.delta1)*1.2])
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
ylim([0, max(opti_coefficients.sigma1)*1.2])
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
ylim([0, max(opti_coefficients.tau2)*1.2])
set(gca, 'TickLabelInterpreter', 'Latex')


% Coefficients lambda
figure(14)
plot(new_time, opti_coefficients.lambda, LineWidth=1.5)
ylabel('Coefficients Values','Interpreter','latex')
title('\textbf{Coefficient $\lambda$}','Interpreter','latex')
grid on
legend('$\lambda$', 'Interpreter','latex', 'Location','northeast')
xlim([new_time(1), new_time(end)])
ylim([0, max(opti_coefficients.lambda)*1.2])
set(gca, 'TickLabelInterpreter', 'Latex')


% Coefficients phi
figure(14)
plot(new_time, opti_coefficients.phi, LineWidth=1.5)
ylabel('Coefficients Values','Interpreter','latex')
title('\textbf{Coefficient $\phi$}','Interpreter','latex')
grid on
legend('$\phi$', 'Interpreter','latex', 'Location','northeast')
xlim([new_time(1), new_time(end)])
ylim([0, max(opti_coefficients.phi)*1.2])
set(gca, 'TickLabelInterpreter', 'Latex')

% Coefficient kappa
figure(15)
plot(new_time, opti_coefficients.kappa, LineWidth=1.5)
ylabel('Coefficients Values','Interpreter','latex')
title('\textbf{Coefficient $\kappa$}','Interpreter','latex')
grid on
legend('$\kappa$', 'Interpreter','latex', 'Location','northeast')
xlim([new_time(1), new_time(end)])
ylim([0, max(opti_coefficients.kappa)*1.2])
set(gca, 'TickLabelInterpreter', 'Latex')

%% Additional Intersting plots,

customColors2 = {   [1, 0.6, 0.2],
                    [1, 0.3, 0.05],
                    [1, 0.2, 0.05],
                    [0.8, 0.2, 0.1],
                    [1, 0.2, 0.05],            
                    [0.8, 0.2, 0.1],
                    [1, 0.3, 0.05],
                    [1, 0.6, 0.2]
                };

for ii = 1:length(policy_dates)-1

    area.x(ii, :) = [policy_dates(ii) policy_dates(ii) policy_dates(ii+1) policy_dates(ii+1)];
    area.y_alpha(ii, :) = [0 max(opti_coefficients.alpha)*1.05 max(opti_coefficients.alpha)*1.05 0];
    area.y_beta(ii, :) = [0 max(opti_coefficients.beta)*1.05 max(opti_coefficients.beta)*1.05 0];
    area.y_gamma(ii, :) = [0 max(opti_coefficients.gamma)*1.05 max(opti_coefficients.gamma)*1.05 0];
    
    area.y_S(ii, :) = [0 max(SIDTTHE_trends.S)*1.05 max(SIDTTHE_trends.S)*1.05 0];
    area.y_I(ii, :) = [0 max(SIDTTHE_trends.I)*1.05 max(SIDTTHE_trends.I)*1.05 0];
    area.y_D(ii, :) = [0 max(SIDTTHE_trends.D)*1.05 max(SIDTTHE_trends.D)*1.05 0];
    area.y_T1(ii, :) = [0 max(DataArray{4,:})*1.05 max(DataArray{4,:})*1.05 0];
    area.y_T2(ii, :) = [0 max(SIDTTHE_trends.T2)*1.05 max(SIDTTHE_trends.T2)*1.05 0];
    area.y_H(ii, :) = [0 max(SIDTTHE_trends.H)*1.05 max(SIDTTHE_trends.H)*1.05 0];
    area.y_E(ii, :) = [0 max(DataArray{7,:})*1.05 max(DataArray{7,:})*1.05 0];
    area.y_V(ii, :) = [0 max(SIDTTHE_trends.V)*1.05 max(SIDTTHE_trends.V)*1.05 0];

end

% Figure of the alpha trend related to policy In italy
figure()
for ii = 1:length(policy_dates)-1
    fill(area.x(ii, :) ,area.y_alpha(1, :), customColors2{ii,1} ,'FaceAlpha',.5,'EdgeColor', 'none','HandleVisibility', 'off')
    hold on
    xline(policy_dates(ii),":",'HandleVisibility', 'off')
    hold on 
end
hold on
plot(new_time, opti_coefficients.alpha,'k','LineWidth',1.5, 'DisplayName', '$\alpha$')
ylabel('Coefficients Values','Interpreter','latex')
title('\textbf{$\alpha$ coefficient}','Interpreter','latex')
grid on
legend('Interpreter','latex','location','southeast')
ylim([0, max(opti_coefficients.alpha)*1.05])
xlim([new_time(1), new_time(end)])
set(gca, 'TickLabelInterpreter', 'Latex')



% Figure of the beta trend related to policy In italy
figure()
for ii = 1:length(policy_dates)-1
    fill(area.x(ii, :) ,area.y_beta(1, :), customColors2{ii,1} ,'FaceAlpha',.5,'EdgeColor', 'none','HandleVisibility', 'off')
    hold on
    xline(policy_dates(ii),":",'HandleVisibility', 'off')
    hold on 
end
hold on
plot(new_time, opti_coefficients.beta,'k','LineWidth',1.5, 'DisplayName', '$\beta$')
ylabel('Coefficients Values','Interpreter','latex')
title('\textbf{$\beta$ coefficient}','Interpreter','latex')
grid on
legend('Interpreter','latex')
xlim([new_time(1), new_time(end)])
ylim([0, max(opti_coefficients.beta)*1.05])
set(gca, 'TickLabelInterpreter', 'Latex')




% Figure of the Gamma trend related to policy In italy
figure()
for ii = 1:length(policy_dates)-1
    fill(area.x(ii, :) ,area.y_gamma(1, :), customColors2{ii,1} ,'FaceAlpha',.5,'EdgeColor', 'none','HandleVisibility', 'off')
    hold on
    xline(policy_dates(ii),":",'HandleVisibility', 'off')
    hold on 
end
hold on
plot(new_time, opti_coefficients.gamma,'k','LineWidth',1.5, 'DisplayName', '$\gamma$')
hold on 
ylabel('Coefficients Values','Interpreter','latex')
title('\textbf{$\gamma$ coefficient}','Interpreter','latex')
grid on
legend('Interpreter','latex','location','southeast')
xlim([new_time(1), new_time(end)])
ylim([0, max(opti_coefficients.gamma)*1.05])
set(gca, 'TickLabelInterpreter', 'Latex')


figure()
% first data set on the primary y-axis
yyaxis left;
plot(new_time, opti_coefficients.gamma,'LineWidth', 1.5);
ylabel('$\gamma$ Coefficients Values','Interpreter','latex')
ax1 = gca;
ax1.YColor = 'k';
% second data set on the secondary y-axis
yyaxis right;
plot(new_time, tests.data/Npop,'-.','LineWidth',1.5,'DisplayName', 'Testing Activity')
ylabel('$\char"0023$ of Tests','Interpreter','latex')
ax2 = gca;
ax2.YColor = 'k';
title('\textbf{Comparison Testing Activity and $\gamma$ Coefficient}','Interpreter','latex')
legend(' $\gamma$', 'Tests','Interpreter','latex')
grid on;
xlim([new_time(1), new_time(end)])

%% Data and fitting in with NPIs
matlabBlue = [0, 0, 0];

% S - data + NPIs
figure(1)
for ii = 1:length(policy_dates)-1
    fill(area.x(ii, :) ,area.y_S(1, :), customColors2{ii,1} ,'FaceAlpha',.5,'EdgeColor', 'none','HandleVisibility', 'off')
    hold on
    xline(policy_dates(ii),":",'HandleVisibility', 'off')
    hold on 
end
hold on
plot(new_time, DataArray{1,:}, LineWidth =1.5, Color=matlabBlue)
ylabel('$\char"0023$ of cases','Interpreter','latex')
title('\textbf{Data - Susceptible}','Interpreter','latex')
legend('Real Data', 'ODE Simulated Data','Interpreter','latex', 'Location','northeast')
xlim([new_time(1), new_time(end)])
ylim([0, max(SIDTTHE_trends.S)*1.05])
set(gca, 'TickLabelInterpreter', 'Latex')


% I - data + NPIs
figure(2)
for ii = 1:length(policy_dates)-1
    fill(area.x(ii, :) ,area.y_I(1, :), customColors2{ii,1} ,'FaceAlpha',.5,'EdgeColor', 'none','HandleVisibility', 'off')
    hold on
    xline(policy_dates(ii),":",'HandleVisibility', 'off')
    hold on 
end
hold on
plot(new_time, DataArray{2,:}, LineWidth =1.5, Color=matlabBlue)
ylabel('$\char"0023$ of cases','Interpreter','latex')
title('\textbf{Data - Infected}','Interpreter','latex')
legend('Real Data', 'ODE Simulated Data','Interpreter','latex', 'Location','northeast')
xlim([new_time(1), new_time(end)])
ylim([0, max(SIDTTHE_trends.I)*1.05])
set(gca, 'TickLabelInterpreter', 'Latex')

% D - data + NPIs
figure(3)
for ii = 1:length(policy_dates)-1
    fill(area.x(ii, :) ,area.y_D(1, :), customColors2{ii,1} ,'FaceAlpha',.5,'EdgeColor', 'none','HandleVisibility', 'off')
    hold on
    xline(policy_dates(ii),":",'HandleVisibility', 'off')
    hold on 
end
hold on
plot(new_time, DataArray{3,:}, LineWidth =1.5, Color=matlabBlue)
ylabel('$\char"0023$ of cases','Interpreter','latex')
title('\textbf{Data - Diagnosed}','Interpreter','latex')
legend('Real Data', 'ODE Simulated Data','Interpreter','latex', 'Location','northeast')
xlim([new_time(1), new_time(end)])
ylim([0, max(SIDTTHE_trends.D)*1.05])
set(gca, 'TickLabelInterpreter', 'Latex')

% T1 - data + NPIs
figure(4)
for ii = 1:length(policy_dates)-1
    fill(area.x(ii, :) ,area.y_T1(1, :), customColors2{ii,1} ,'FaceAlpha',.5,'EdgeColor', 'none','HandleVisibility', 'off')
    hold on
    xline(policy_dates(ii),":",'HandleVisibility', 'off')
    hold on 
end
hold on
plot(new_time, DataArray{4,:}, LineWidth =1.5, Color=matlabBlue)
ylabel('$\char"0023$ of cases','Interpreter','latex')
title('\textbf{Data - Hospitalised}','Interpreter','latex')
legend('Real Data', 'ODE Simulated Data','Interpreter','latex', 'Location','northeast')
xlim([new_time(1), new_time(end)])
ylim([0, max(DataArray{4,:})*1.05])
set(gca, 'TickLabelInterpreter', 'Latex')

% T2 - data + NPIs
figure(5)
for ii = 1:length(policy_dates)-1
    fill(area.x(ii, :) ,area.y_T2(1, :), customColors2{ii,1} ,'FaceAlpha',.5,'EdgeColor', 'none','HandleVisibility', 'off')
    hold on
    xline(policy_dates(ii),":",'HandleVisibility', 'off')
    hold on 
end
hold on
plot(new_time, DataArray{5,:}, LineWidth =1.5, Color=matlabBlue)
ylabel('$\char"0023$ of cases','Interpreter','latex')
title('\textbf{Data - ICUs}','Interpreter','latex')
legend('Real Data', 'ODE Simulated Data','Interpreter','latex', 'Location','northeast')
xlim([new_time(1), new_time(end)])
ylim([0, max(SIDTTHE_trends.T2)*1.05])
set(gca, 'TickLabelInterpreter', 'Latex')

% H - data + NPIs
figure(6)
for ii = 1:length(policy_dates)-1
    fill(area.x(ii, :) ,area.y_H(1, :), customColors2{ii,1} ,'FaceAlpha',.5,'EdgeColor', 'none','HandleVisibility', 'off')
    hold on
    xline(policy_dates(ii),":",'HandleVisibility', 'off')
    hold on 
end
hold on
plot(new_time, DataArray{6,:}, LineWidth =1.5, Color=matlabBlue)
ylabel('$\char"0023$ of cases','Interpreter','latex')
title('\textbf{Data - Healed}','Interpreter','latex')
legend('Real Data', 'ODE Simulated Data','Interpreter','latex', 'Location','southeast')
xlim([new_time(1), new_time(end)])
ylim([0, max(SIDTTHE_trends.H)*1.05])
set(gca, 'TickLabelInterpreter', 'Latex')

% E - data + NPIs
figure(7)
for ii = 1:length(policy_dates)-1
    fill(area.x(ii, :) ,area.y_E(1, :), customColors2{ii,1} ,'FaceAlpha',.5,'EdgeColor', 'none','HandleVisibility', 'off')
    hold on
    xline(policy_dates(ii),":",'HandleVisibility', 'off')
    hold on 
end
hold on
plot(new_time, DataArray{7,:}, LineWidth =1.5, Color=matlabBlue)
ylabel('$\char"0023$ of cases','Interpreter','latex')
title('\textbf{Data - Deceased}','Interpreter','latex')
legend('Real Data', 'ODE Simulated Data','Interpreter','latex', 'Location','southeast')
xlim([new_time(1), new_time(end)])
ylim([0, max(DataArray{7,:})*1.05])
set(gca, 'TickLabelInterpreter', 'Latex')

% V - data + NPIs
figure(8)
for ii = 1:length(policy_dates)-1
    fill(area.x(ii, :) ,area.y_V(1, :), customColors2{ii,1} ,'FaceAlpha',.5,'EdgeColor', 'none','HandleVisibility', 'off')
    hold on
    xline(policy_dates(ii),":",'HandleVisibility', 'off')
    hold on 
end
hold on
plot(new_time, DataArray{8,:}, LineWidth =1.5, Color=matlabBlue)
ylabel('$\char"0023$ of cases','Interpreter','latex')
title('\textbf{Data - Vaccinated}','Interpreter','latex')
legend('Real Data', 'ODE Simulated Data','Interpreter','latex', 'Location','southeast')
xlim([new_time(1), new_time(end)])
ylim([0, max(SIDTTHE_trends.V)*1.05])
set(gca, 'TickLabelInterpreter', 'Latex')


%% Variants and coefficients 

load("/Users/marcodelloro/Desktop/Thesis/MSc-Thesis-TUe/Data_Collection/Smooth_Variants.mat")

% Alpha and Variants
figure()
yyaxis right;
fill([SmoothVar.date; flip(SmoothVar.date)], [SmoothVar.SC2 zeros(1,57)], [0 0.4470 0.7410], 'LineStyle', '-','FaceAlpha',.1,'EdgeColor', [0 0.4470 0.7410])
hold on
fill([SmoothVar.date; flip(SmoothVar.date)], [SmoothVar.alpha zeros(1,57)], [0.8500 0.3250 0.0980], 'LineStyle', '-', 'FaceAlpha',.1,'EdgeColor', [0.8500 0.3250 0.0980])
hold on 
fill([SmoothVar.date; flip(SmoothVar.date)], [SmoothVar.gamma zeros(1,57)], [0.9290 0.6940 0.1250], 'LineStyle', '-', 'FaceAlpha',.1,'EdgeColor', [0.9290 0.6940 0.1250])
hold on
fill([SmoothVar.date; flip(SmoothVar.date)], [SmoothVar.delta zeros(1,57)], [0.4940 0.1840 0.5560], 'LineStyle', '-', 'FaceAlpha',.1,'EdgeColor', [0.4940 0.1840 0.5560])
ax2 = gca;
ax2.YColor = 'k';
yyaxis left;
plot(new_time, opti_coefficients.alpha, 'k','LineWidth', 2);
ylabel('$\beta$ Coefficient Values','Interpreter','latex')
ax1 = gca;
ax1.YColor = 'k';
title('\textbf{}','Interpreter','latex')
title('\textbf{SARS-CoV-2 Variants - Coefficient $\alpha$}','Interpreter','latex')
legend('$\alpha$','SARS-CoV-2', 'B.1.1.7 - ALPHA VARIANT', 'P.1 - GAMMA VARIANT', ' B.1.617.2 - DELTA VARIANT', 'Interpreter','latex','Location','northeastoutside')
grid on
xlim([SmoothVar.date(1), SmoothVar.date(end)])
ylim([0, max(opti_coefficients.alpha)*1.2])

% Beta and Variants
figure()
yyaxis right;
fill([SmoothVar.date; flip(SmoothVar.date)], [SmoothVar.SC2 zeros(1,57)], [0 0.4470 0.7410], 'LineStyle', '-','FaceAlpha',.1,'EdgeColor', [0 0.4470 0.7410])
hold on
fill([SmoothVar.date; flip(SmoothVar.date)], [SmoothVar.alpha zeros(1,57)], [0.8500 0.3250 0.0980], 'LineStyle', '-', 'FaceAlpha',.1,'EdgeColor', [0.8500 0.3250 0.0980])
hold on 
fill([SmoothVar.date; flip(SmoothVar.date)], [SmoothVar.gamma zeros(1,57)], [0.9290 0.6940 0.1250], 'LineStyle', '-', 'FaceAlpha',.1,'EdgeColor', [0.9290 0.6940 0.1250])
hold on
fill([SmoothVar.date; flip(SmoothVar.date)], [SmoothVar.delta zeros(1,57)], [0.4940 0.1840 0.5560], 'LineStyle', '-', 'FaceAlpha',.1,'EdgeColor', [0.4940 0.1840 0.5560])
ax2 = gca;
ax2.YColor = 'k';
yyaxis left;
plot(new_time, opti_coefficients.beta, 'k','LineWidth', 2);
ylabel('$\beta$ Coefficient Values','Interpreter','latex')
ax1 = gca;
ax1.YColor = 'k';
title('\textbf{}','Interpreter','latex')
title('\textbf{SARS-CoV-2 Variants - Coefficient $\beta$}','Interpreter','latex')
legend('$\beta$','SARS-CoV-2', 'B.1.1.7 - ALPHA VARIANT', 'P.1 - GAMMA VARIANT', ' B.1.617.2 - DELTA VARIANT', 'Interpreter','latex','Location','northeastoutside')
grid on
xlim([SmoothVar.date(1), SmoothVar.date(end)])
ylim([0, max(opti_coefficients.beta)*1.2])


% Deltas and Variants
figure()
yyaxis right;
fill([SmoothVar.date; flip(SmoothVar.date)], [SmoothVar.SC2 zeros(1,57)], [0 0.4470 0.7410], 'LineStyle', '-','FaceAlpha',.1,'EdgeColor', [0 0.4470 0.7410])
hold on
fill([SmoothVar.date; flip(SmoothVar.date)], [SmoothVar.alpha zeros(1,57)], [0.8500 0.3250 0.0980], 'LineStyle', '-', 'FaceAlpha',.1,'EdgeColor', [0.8500 0.3250 0.0980])
hold on 
fill([SmoothVar.date; flip(SmoothVar.date)], [SmoothVar.gamma zeros(1,57)], [0.9290 0.6940 0.1250], 'LineStyle', '-', 'FaceAlpha',.1,'EdgeColor', [0.9290 0.6940 0.1250])
hold on
fill([SmoothVar.date; flip(SmoothVar.date)], [SmoothVar.delta zeros(1,57)], [0.4940 0.1840 0.5560], 'LineStyle', '-', 'FaceAlpha',.1,'EdgeColor', [0.4940 0.1840 0.5560])
ax2 = gca;
ax2.YColor = 'k';
yyaxis left;
plot(new_time, opti_coefficients.delta1, 'k-', LineWidth=2)
hold on
plot(new_time, opti_coefficients.delta2, 'r-', LineWidth=2)
ylabel('$\delta$ Coefficient Values','Interpreter','latex')
ax1 = gca;
ax1.YColor = 'k';
title('\textbf{}','Interpreter','latex')
title('\textbf{SARS-CoV-2 Variants - Coefficient $\delta$}','Interpreter','latex')
legend('$\delta_1$','$\delta_2$','SARS-CoV-2', 'B.1.1.7 - ALPHA VARIANT', 'P.1 - GAMMA VARIANT', ' B.1.617.2 - DELTA VARIANT', 'Interpreter','latex','Location','northeastoutside')
grid on
xlim([SmoothVar.date(1), SmoothVar.date(end)])
ylim([0, max(opti_coefficients.delta1)*1.2])

% Sigmas and Variants
figure()
yyaxis right;
fill([SmoothVar.date; flip(SmoothVar.date)], [SmoothVar.SC2 zeros(1,57)], [0 0.4470 0.7410], 'LineStyle', '-','FaceAlpha',.1,'EdgeColor', [0 0.4470 0.7410])
hold on
fill([SmoothVar.date; flip(SmoothVar.date)], [SmoothVar.alpha zeros(1,57)], [0.8500 0.3250 0.0980], 'LineStyle', '-', 'FaceAlpha',.1,'EdgeColor', [0.8500 0.3250 0.0980])
hold on 
fill([SmoothVar.date; flip(SmoothVar.date)], [SmoothVar.gamma zeros(1,57)], [0.9290 0.6940 0.1250], 'LineStyle', '-', 'FaceAlpha',.1,'EdgeColor', [0.9290 0.6940 0.1250])
hold on
fill([SmoothVar.date; flip(SmoothVar.date)], [SmoothVar.delta zeros(1,57)], [0.4940 0.1840 0.5560], 'LineStyle', '-', 'FaceAlpha',.1,'EdgeColor', [0.4940 0.1840 0.5560])
ax2 = gca;
ax2.YColor = 'k';
yyaxis left;
plot(new_time, opti_coefficients.sigma1, 'k-', LineWidth=2)
hold on
plot(new_time, opti_coefficients.sigma2, 'r-', LineWidth=2)
ylabel('$\sigma$ Coefficient Values','Interpreter','latex')
ax1 = gca;
ax1.YColor = 'k';
title('\textbf{}','Interpreter','latex')
title('\textbf{SARS-CoV-2 Variants - Coefficient $\sigma$}','Interpreter','latex')
legend('$\sigma_1$','$\sigma_2$','SARS-CoV-2', 'B.1.1.7 - ALPHA VARIANT', 'P.1 - GAMMA VARIANT', ' B.1.617.2 - DELTA VARIANT', 'Interpreter','latex','Location','northeastoutside')
grid on
xlim([SmoothVar.date(1), SmoothVar.date(end)])
ylim([0, max(opti_coefficients.sigma1)*1.2])


% Taus and Variants
figure()
yyaxis right;
fill([SmoothVar.date; flip(SmoothVar.date)], [SmoothVar.SC2 zeros(1,57)], [0 0.4470 0.7410], 'LineStyle', '-','FaceAlpha',.1,'EdgeColor', [0 0.4470 0.7410])
hold on
fill([SmoothVar.date; flip(SmoothVar.date)], [SmoothVar.alpha zeros(1,57)], [0.8500 0.3250 0.0980], 'LineStyle', '-', 'FaceAlpha',.1,'EdgeColor', [0.8500 0.3250 0.0980])
hold on 
fill([SmoothVar.date; flip(SmoothVar.date)], [SmoothVar.gamma zeros(1,57)], [0.9290 0.6940 0.1250], 'LineStyle', '-', 'FaceAlpha',.1,'EdgeColor', [0.9290 0.6940 0.1250])
hold on
fill([SmoothVar.date; flip(SmoothVar.date)], [SmoothVar.delta zeros(1,57)], [0.4940 0.1840 0.5560], 'LineStyle', '-', 'FaceAlpha',.1,'EdgeColor', [0.4940 0.1840 0.5560])
ax2 = gca;
ax2.YColor = 'k';
yyaxis left;
plot(new_time, opti_coefficients.tau1, 'k-', LineWidth=2)
hold on
plot(new_time, opti_coefficients.tau2, 'r-', LineWidth=2)
ylabel('$\tau$ Coefficient Values','Interpreter','latex')
ax1 = gca;
ax1.YColor = 'k';
title('\textbf{}','Interpreter','latex')
title('\textbf{SARS-CoV-2 Variants - Coefficient $\tau$}','Interpreter','latex')
legend('$\tau_1$','$\tau_2$','SARS-CoV-2', 'B.1.1.7 - ALPHA VARIANT', 'P.1 - GAMMA VARIANT', ' B.1.617.2 - DELTA VARIANT', 'Interpreter','latex','Location','northeastoutside')
grid on
xlim([SmoothVar.date(1), SmoothVar.date(end)])
ylim([0, max(opti_coefficients.tau2)*1.2])


% Lambda and Variants
figure()
yyaxis right;
fill([SmoothVar.date; flip(SmoothVar.date)], [SmoothVar.SC2 zeros(1,57)], [0 0.4470 0.7410], 'LineStyle', '-','FaceAlpha',.1,'EdgeColor', [0 0.4470 0.7410])
hold on
fill([SmoothVar.date; flip(SmoothVar.date)], [SmoothVar.alpha zeros(1,57)], [0.8500 0.3250 0.0980], 'LineStyle', '-', 'FaceAlpha',.1,'EdgeColor', [0.8500 0.3250 0.0980])
hold on 
fill([SmoothVar.date; flip(SmoothVar.date)], [SmoothVar.gamma zeros(1,57)], [0.9290 0.6940 0.1250], 'LineStyle', '-', 'FaceAlpha',.1,'EdgeColor', [0.9290 0.6940 0.1250])
hold on
fill([SmoothVar.date; flip(SmoothVar.date)], [SmoothVar.delta zeros(1,57)], [0.4940 0.1840 0.5560], 'LineStyle', '-', 'FaceAlpha',.1,'EdgeColor', [0.4940 0.1840 0.5560])
ax2 = gca;
ax2.YColor = 'k';
yyaxis left;
plot(new_time, opti_coefficients.lambda, 'k','LineWidth', 2);
ylabel('$\lambda$ Coefficient Values','Interpreter','latex')
ax1 = gca;
ax1.YColor = 'k';
title('\textbf{}','Interpreter','latex')
title('\textbf{SARS-CoV-2 Variants - Coefficient $\lambda$}','Interpreter','latex')
legend('$\lambda$','SARS-CoV-2', 'B.1.1.7 - ALPHA VARIANT', 'P.1 - GAMMA VARIANT', ' B.1.617.2 - DELTA VARIANT', 'Interpreter','latex','Location','northeastoutside')
grid on
xlim([SmoothVar.date(1), SmoothVar.date(end)])
ylim([0, max(opti_coefficients.lambda)*1.2])

% Epsi and Variants
figure()
yyaxis right;
fill([SmoothVar.date; flip(SmoothVar.date)], [SmoothVar.SC2 zeros(1,57)], [0 0.4470 0.7410], 'LineStyle', '-','FaceAlpha',.1,'EdgeColor', [0 0.4470 0.7410])
hold on
fill([SmoothVar.date; flip(SmoothVar.date)], [SmoothVar.alpha zeros(1,57)], [0.8500 0.3250 0.0980], 'LineStyle', '-', 'FaceAlpha',.1,'EdgeColor', [0.8500 0.3250 0.0980])
hold on 
fill([SmoothVar.date; flip(SmoothVar.date)], [SmoothVar.gamma zeros(1,57)], [0.9290 0.6940 0.1250], 'LineStyle', '-', 'FaceAlpha',.1,'EdgeColor', [0.9290 0.6940 0.1250])
hold on
fill([SmoothVar.date; flip(SmoothVar.date)], [SmoothVar.delta zeros(1,57)], [0.4940 0.1840 0.5560], 'LineStyle', '-', 'FaceAlpha',.1,'EdgeColor', [0.4940 0.1840 0.5560])
ax2 = gca;
ax2.YColor = 'k';
yyaxis left;
plot(new_time, opti_coefficients.epsi, 'k','LineWidth', 2);
ylabel('$\epsilon$ Coefficient Values','Interpreter','latex')
ax1 = gca;
ax1.YColor = 'k';
title('\textbf{}','Interpreter','latex')
title('\textbf{SARS-CoV-2 Variants - Coefficient $\epsilon$}','Interpreter','latex')
legend('$\epsilon$','SARS-CoV-2', 'B.1.1.7 - ALPHA VARIANT', 'P.1 - GAMMA VARIANT', ' B.1.617.2 - DELTA VARIANT', 'Interpreter','latex','Location','northeastoutside')
grid on
xlim([SmoothVar.date(1), SmoothVar.date(end)])