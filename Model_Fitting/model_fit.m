%% Version 2 - No Epsilon 1 value, no transmission between T1 and T2
% In this version model_fit we will keep the main model 

% This version uses a smoothing of the coeffcients by minimizing the
% difference between the coefficients - best model yet

clc
clear all 

load('/Users/marcodelloro/Desktop/Thesis/MSc-Thesis-TUe/Data_Colection/SIDTTHE_data.mat');
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
H_data = SIDTTHE_data{6,1}.data / Npop;     % T2 group - Healed Population **
E_data = SIDTTHE_data{5,1}.data / Npop;     % E group - Expired (Deceased) Population

S_data = ones(length(I_data),1) - (I_data + D_data + T1_data + T2_data + H_data + E_data); 

% ** the groups I and H are based on strong assumptions and dark numbers estimates

DataArray = {S_data, I_data, D_data, T1_data, T2_data, H_data, E_data};

N = 399;

%% Choice of the sampling method 

% Create a dialog box with three buttons
dlgTitle = 'Sampling Method Selection';
prompt = 'Select a sampling method:';
options = {'57 points sampling', 'Latin Hypercube sampling', 'Full sampling'};
defaultOption = options{1};

selectedOption = questdlg(prompt, dlgTitle, options{:}, defaultOption);

% Check which option the user selected
if isempty(selectedOption)
    % User canceled the dialog
    disp('Sampling method selection canceled.');
else
    disp(['Selected option: ' selectedOption]);


    % Perform the action based on the selected option - 57 qweekly sampling action
    if strcmp(selectedOption, '57 points sampling')
        
        data = horzcat(DataArray{:})';

    elseif strcmp(selectedOption, 'Latin Hypercube sampling')
        
        perc = 0.5;  % percentage of points we want to sample from the initial dataset
        sampled_points = round(N*perc);
        disp(['Sampling method LH - ' num2str(sampled_points) ' points.']);
        
        % New time vector (daily data)
        new_time = linspace(SIDTTHE_data{1, 1}.date(1), SIDTTHE_data{1, 1}.date(end), N); % 399 data points (daily data)
        
        % Iterate over each cell in DataArray
        for ii = 1:numel(DataArray)
        
            original_data = DataArray{1, ii};
        
            upsampled_DataArray{ii,1} = interp1(SIDTTHE_data{1, 1}.date, original_data, new_time);
        
        end
        
        N_new = sort( round(lhsdesign(sampled_points, 1) * N) );
        
        for ii=1:7
            originalData = upsampled_DataArray{ii};
            lhs_DataArray{ii,:} = originalData(N_new);
        end

        data = vertcat(lhs_DataArray{:});

    elseif strcmp(selectedOption, 'Full sampling')
  
        % Upsampling of the dataset using Interp1, so that we can handle better the divided day by day 
        
        N = 399; % if we are using the DAILY time discretization, then N = 399
                 % else is N = 57 (weekly)
        
        % New time vector (daily data)
        new_time = linspace(SIDTTHE_data{1, 1}.date(1), SIDTTHE_data{1, 1}.date(end), N); % 399 data points (daily data)
        
        % Iterate over each cell in DataArray
        for ii = 1:numel(DataArray)
        
            original_data = DataArray{1, ii};
        
            upsampled_DataArray{ii,1} = interp1(SIDTTHE_data{1, 1}.date, original_data, new_time);
        
        end

         data = vertcat(upsampled_DataArray{:});

    else
        disp('Invalid selection.');
    end
end


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

opti.set_initial(S,ones(N,1) * S_data(1));
opti.set_initial(I,ones(N,1) * I_data(1));
opti.set_initial(D,ones(N,1) * D_data(1));
opti.set_initial(T1,ones(N,1)* T1_data(1));
opti.set_initial(T2,ones(N,1)* T2_data(1));
opti.set_initial(H,ones(N,1) * H_data(1));
opti.set_initial(E,ones(N,1) * E_data(1));

opti.set_initial(alpha,ones(N,1) * 0.513);          % All coeff. are obtained from the mean of the values provided by G.Giordano in original SIDARTHE paper
opti.set_initial(beta,ones(N,1) * 0.011); 
opti.set_initial(gamma,ones(N,1) * 0.2405);
opti.set_initial(delta1,ones(N,1) * 0.01);
opti.set_initial(delta2,ones(N,1) * 0.008);

opti.set_initial(sigma1,ones(N,1) * 0.03);
opti.set_initial(sigma2,ones(N,1) * 0.065);
opti.set_initial(tau1,ones(N,1) * 0.015);
opti.set_initial(tau2,ones(N,1) * 0.005);
opti.set_initial(lambda,ones(N,1) * 0.0596);

%% Constraints for the parameters

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
opti.subject_to(delta1(1:end) <= 0.02);

opti.subject_to(delta2(1:end) >= 0);     % bound on delta2 value
opti.subject_to(delta2(1:end) <= 0.001);

% epsi parameter - bound on epsilon value (different bound for different months)

opti.subject_to(epsi(1:end) >= 0);    % bound on epsi value
opti.subject_to(epsi(1:end) <= 0.007);

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

% Initial Conditions constraints

opti.subject_to(S(1,1) == data(1,1));
opti.subject_to(I(1,1) == data(2,1));
opti.subject_to(D(1,1) == data(3,1));
opti.subject_to(T1(1,1) == data(4,1));
opti.subject_to(T2(1,1) == data(5,1));
opti.subject_to(H(1,1) == data(6,1));
opti.subject_to(E(1,1) == data(7,1));

opti.subject_to(S + I + D + T1 + T2 + H + E == 1);

opti.subject_to(S >= 0);
opti.subject_to(I >= 0);
opti.subject_to(D >= 0);
opti.subject_to(T1 >= 0);
opti.subject_to(T2 >= 0);
opti.subject_to(H >= 0);
opti.subject_to(E >= 0);

opti.subject_to(alpha >= 2.5*beta); % bound on contagion parameter between alpha and beta

% Constraints based on assumptions and previous Knowledge                                                                                        

opti.subject_to(tau1 >= tau2);
opti.subject_to(sigma2 >= sigma1);
opti.subject_to(delta1 > 2 * delta2);


%% RK45 simulation of the differential equations

% Specify system dynamics - creation of the handle functions
% Version 2
f = @(X, coefs)   [           -X(1)*(coefs(1)*X(2) + coefs(2)*X(3));                                                 % dS/dt = S_dot 
                              ( X(1)*(coefs(1)*X(2) + coefs(2)*X(3)) ) - ( (coefs(3) + coefs(11))*X(2) );             % dI/dt = I_dot
                              (X(2)*coefs(3)) - ( X(3)*(coefs(11) + coefs(4) + coefs(5)) );                           % dD/dt = D_dot
                              ( coefs(4)*X(3) ) - ( (coefs(7) + coefs(9) + coefs(6))*X(4) );                            % dT1/dt = T1_dot
                              ( coefs(5)*X(3) ) - ( (coefs(8) + coefs(10))*X(5) ) + ( coefs(6)*X(4) );                     % dT2/dt = T2_dot
                              ( (X(2) + X(3))*coefs(11) ) + ( X(4)*coefs(7) ) + ( X(5)*coefs(8) ) ;                       % dH/dt = H_dot
                              ( X(4)*coefs(9) ) + ( X(5)*coefs(10) ) ];  

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

% 8 different areas will be covered
policy_idx = [1 40 69 116 141 190 242 294 399];   % Values taken from italian policies applied for the pandemic

for ii=1:length(policy_idx)
    policy_dates(ii) = [new_time(policy_idx(ii))];
end

data_obj = horzcat(data(:));
X_obj = [];

% differentiation of X_obj depending on the different sampling methods
if strcmp(selectedOption, '57 points sampling')
    idx = 7:7:N;
    idx = [1 idx];

    for ii=1:7
        X_obj = [X_obj horzcat(X(ii,idx))]';
    end

elseif strcmp(selectedOption, 'Latin Hypercube sampling')
    
    for ii=1:7
        X_obj = [X_obj horzcat(X(ii,N_new))];
    end
    
    X_obj = X_obj';

elseif strcmp(selectedOption, 'Full sampling')
    
    X_obj = horzcat(X(:));

end

% building the actual COST FUNCTIONS FOR THE CONSTRAINT

coefs_matrix_obj = [    sum( ( (diff(coefs(1,:)))./0.6).^2 );    
                        sum( ( (diff(coefs(2,:)))./0.015).^2 );
                        sum( ( (diff(coefs(3,:)))./0.5).^2 );
                        sum( ( (diff(coefs(4,:)))./0.02).^2 );
                        sum( ( (diff(coefs(5,:)))./0.03).^2 );
                        sum( ( (diff(coefs(6,:)))./0.6).^2 );
                        sum( ( (diff(coefs(7,:)))./0.03).^2 );
                        sum( ( (diff(coefs(8,:)))./0.03).^2 );
                        sum( ( (diff(coefs(9,:)))./0.03).^2 );
                        sum( ( (diff(coefs(10,:)))./0.05).^2 );
                        sum( ( (diff(coefs(11,:)))./0.5 ).^2)     ];

cost_matrix_obj =  coefs_matrix_obj(1) + coefs_matrix_obj(2) + coefs_matrix_obj(3) + coefs_matrix_obj(4) + coefs_matrix_obj(5) +...
                   coefs_matrix_obj(6) + coefs_matrix_obj(7) +  coefs_matrix_obj(8) + coefs_matrix_obj(9) + coefs_matrix_obj(10) + coefs_matrix_obj(11) ;

% Cumulative cost between the data (INTEGRAL COST)

integral_cost1= 0; % initialization of the integral cost for Hosipitalised
integral_cost2= 0; % initialization of the integral cost for ICUs

for ii = 1:N
    step_cost1 = (data(4,ii) - X(4,ii)).^2;
    step_cost2 = (data(5,ii) - X(5,ii)).^2;
    integral_cost1 = integral_cost1 + step_cost1;
    integral_cost2 = integral_cost2 + step_cost2;
end

% Definitions of some weights for the optimization

w1 = 5;     % weight for the difference between data and model
w2 = 0.05;   % weight for the coefficient smoothness
w3 = 4;     % weight for the difference between data and model (integral cost Hosipitalised)
w4 = 5;   % weight for the difference between data and model (integral cost ICUs)

obj = sum((data_obj - X_obj).^2) * w1  + cost_matrix_obj * w2 + integral_cost1 * w3 + integral_cost2 * w4 ; 

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
close all

% Coefficients alpha and beta
figure(7)
plot(new_time, opti_coefficients.alpha, LineWidth=1.5)
ylabel('Coefficients Values','Interpreter','latex')
title('\textbf{Coefficient $\alpha$}','Interpreter','latex')
grid on
legend('$\beta$','Interpreter','latex', 'Location','northeast')
xlim([new_time(1), new_time(end)])
set(gca, 'TickLabelInterpreter', 'Latex')


% Coefficients alpha and beta
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


% Coefficient epsi
figure(11)
plot(new_time, opti_coefficients.epsi, LineWidth=1.5)
ylabel('Coefficients Values','Interpreter','latex')
title('\textbf{Coefficient $\epsilon_2$}','Interpreter','latex')
grid on
legend('$\epsilon_2$','Interpreter','latex', 'Location','northeast')
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
legend('Interpreter','latex')
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
ylim([4.8e-3, max(opti_coefficients.beta)*1.05])
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
ylabel('Coefficients Values','Interpreter','latex')
title('\textbf{$\gamma$ coefficient}','Interpreter','latex')
grid on
legend('Interpreter','latex')
xlim([new_time(1), new_time(end)])
ylim([0, max(opti_coefficients.gamma)*1.05])
set(gca, 'TickLabelInterpreter', 'Latex')

%% Means of the coefficients 

% Averaging of the data by week
window = 7; % weekly average

for ii = length(opti_coefficients.alpha):-window:1

    avg_coefs.alpha(ii)= mean(opti_coefficients.alpha(ii - window +1:ii));      
    avg_coefs.beta(ii) = mean(opti_coefficients.beta(ii - window +1:ii));     
    avg_coefs.gamma(ii) = mean(opti_coefficients.gamma(ii - window +1:ii));
    avg_coefs.delta1(ii) = mean(opti_coefficients.delta1(ii - window +1:ii));         
    avg_coefs.delta2(ii) = mean(opti_coefficients.delta2(ii - window +1:ii));
    avg_coefs.epsi(ii) = mean(opti_coefficients.epsi(ii - window +1:ii));
    avg_coefs.sigma1(ii) = mean(opti_coefficients.sigma1(ii - window +1:ii));
    avg_coefs.sigma2(ii) = mean(opti_coefficients.sigma2(ii - window +1:ii));
    avg_coefs.tau1(ii) = mean(opti_coefficients.tau1(ii - window +1:ii));
    avg_coefs.tau2(ii) = mean(opti_coefficients.tau2(ii - window +1:ii));
    avg_coefs.lambda(ii) = mean(opti_coefficients.lambda(ii - window +1:ii));

end

avg_coefs.alpha = (avg_coefs.alpha);
avg_coefs.beta = (avg_coefs.beta);
avg_coefs.gamma = (avg_coefs.gamma);
avg_coefs.delta1 = (avg_coefs.delta1);
avg_coefs.delta2 = (avg_coefs.delta2);
avg_coefs.epsi = (avg_coefs.epsi);
avg_coefs.sigma1 = (avg_coefs.sigma1);
avg_coefs.sigma2 = (avg_coefs.sigma2);
avg_coefs.tau1 = (avg_coefs.tau1);
avg_coefs.tau2 = (avg_coefs.tau2);
avg_coefs.lambda = (avg_coefs.lambda);

% create a vector with constant values (only for alpha and beta)

nonZeroAlpha = avg_coefs.alpha(avg_coefs.alpha ~= 0);
nonZeroBeta = avg_coefs.beta(avg_coefs.beta ~= 0);

avg_coefs.alpha = repelem(nonZeroAlpha, 7);
avg_coefs.beta = repelem(nonZeroBeta, 7);


figure()
plot(new_time, opti_coefficients.alpha, LineWidth=1.5)
hold on
plot(new_time, avg_coefs.alpha, LineWidth=1.5)
ylabel('Coefficients Values','Interpreter','latex')
title('\textbf{$\alpha$ coefficient}','Interpreter','latex')
grid on
legend('$\alpha$','$\alpha$ averaged','Interpreter','latex', 'Location','southeast')
xlim([new_time(1), new_time(end)])
set(gca, 'TickLabelInterpreter', 'Latex')

figure()
plot(new_time, opti_coefficients.beta, LineWidth=1.5)
hold on
plot(new_time, avg_coefs.beta, LineWidth=1.5)
ylabel('Coefficients Values','Interpreter','latex')
title('\textbf{$\beta$ coefficient}','Interpreter','latex')
grid on
legend('$\beta$','$\beta$ averaged','Interpreter','latex', 'Location','southeast')
xlim([new_time(1), new_time(end)])
set(gca, 'TickLabelInterpreter', 'Latex')

