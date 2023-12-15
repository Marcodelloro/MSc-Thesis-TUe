%% EKF - estimation of unknown parameters

clear all
close all
clc

load('/Users/marcodelloro/Desktop/Thesis/MSc-Thesis-TUe/Data_Collection/SIDTTHE_data.mat');
load('/Users/marcodelloro/Desktop/Thesis/MSc-Thesis-TUe/Data_Collection/Tests_data.mat');
load('/Users/marcodelloro/Desktop/Thesis/MSc-Thesis-TUe/Model_Fitting/EKF_coefs.mat');

addpath('/Users/marcodelloro/Downloads/casadi-3.6.3-osx64-matlab2018b')
import casadi.*;
opti = casadi.Opti();

date = SIDTTHE_data{1,1}.date;
Npop = 59240329; % Total Population of Italy

I_data = SIDTTHE_data{1,1}.data / Npop;      % I group - Infected Population
D_data = SIDTTHE_data{2,1}.data / Npop;      % D group - Diagnosed Population
T1_data = SIDTTHE_data{3,1}.data / Npop;     % T1 group - Hospitalised (Lightly Threatned) Population
T2_data = SIDTTHE_data{4,1}.data / Npop;     % T2 group - ICUs (Severly Threatned) Population
H_data = SIDTTHE_data{6,1}.data / Npop;      % T2 group - Healed Population
E_data = SIDTTHE_data{5,1}.data  / Npop;     % E group - Expired (Deceased) Population
S_data = ones(length(I_data),1)' - (I_data + D_data + T1_data + T2_data + H_data + E_data );   % S group - Susceptible Population

%% Algorithm initialization
nstates = 17; % number of states - already augmented
nmeasures = 5; % number of measurements - already augmented

% XstateA = fi([], 1, 16) S(alpha,beta) (X(1))
%            I(alpha,beta) (X(2))
%            D             (X(3))
%            T1            (X(4))
%            T2            (X(5))
%            H             (X(6))
%            E             (X(7))
%            alpha         (X(8))
%            beta          (X(9))
%            gamma         (X(10))
%            delta1        (X(11))
%            delta2        (X(12))
%            epsi          (X(13))
%            sigma1        (X(14))
%            sigma2        (X(15))
%            tau1          (X(16))
%            tau2          (X(17))

% First guess estimate of the augmented state
alpha_0 = 0.6;
beta_0 = 0.02;
gamma_0 = 0.5;
delta1_0 = 0.002;
delta2_0 = 5e-4;
epsi_0 = 5e-4;
sigma1_0 = 0.03;
sigma2_0 = 0.01;
tau1_0 = 0.03;
tau2_0 = 0.008;

x0A = [S_data(1); I_data(1); D_data(1); T1_data(1); T2_data(1); H_data(1); E_data(1); alpha_0; beta_0; gamma_0; delta1_0; delta2_0; epsi_0; sigma1_0; sigma2_0; tau1_0; tau2_0 ];

% Initial covariance matrix P0
P0 = diag([50 50 5 5 5 5 5 50 50 50 50 50 50 50 50 50 50]); % assumed high covariance only on the PREDICTED UNMEASURED dynamics

% Noise & Input covariance

n = 0.0001; % Noise amplitude in the measurements - measurement without noise

R = n*eye(nmeasures);  % Artificial disturbances on alpha dynamics

Zmeas = [D_data; T1_data; T2_data; H_data; E_data ];

% Process noise / disturbanceson the dynamics of the coefficents

w_alpha = 0.01; % Artificial disturbances on dynamics of coefficients
w_beta = 0.05;
w_gamma = 0.01;
w_delta1 = 0.05;
w_delta2 = 0.01;
w_epsi = 0.05;
w_sigma1 = 0.01;
w_sigma2= 0.05;
w_tau1 = 0.01;
w_tau2 = 0.05;

w_dist = [w_alpha w_beta w_gamma w_delta1 w_delta2 w_epsi w_sigma1 w_sigma2 w_tau1 w_tau2];

for ii = 1:399

    Zmeas_k = Zmeas(:,ii);

    [Xk_plus, Pk_plus,Ob,z,Abar,Bbar,Cbar,T,k,A,C] = sidttheEKF(x0A,P0,R,w_dist,Zmeas_k);

    x0A = Xk_plus;
    P0 = Pk_plus;

    Xtot(:,ii) = Xk_plus;

end

% Saving data to a .mat file
save('EKF_file.mat', 'Ob', 'z', 'Abar', 'Cbar', 'T', 'k', 'A', 'C');

%% Coefficients and States trend plot

close all

% Susceptible - S
figure(1)
plot(date, S_data, LineWidth=1.5)
hold on
plot(date, Xtot(1,:), LineWidth=1.5)
ylabel('$\char"0023$ of cases','Interpreter','latex')
title('\textbf{Curve Fitting - Susceptible}','Interpreter','latex')
grid on
legend('Real Data', 'EKF Estimated Data','Interpreter','latex', 'Location','northeast')
xlim([date(1), date(end)])
set(gca, 'TickLabelInterpreter', 'Latex')

% Infected - I
figure(2)
plot(date, I_data, LineWidth=1.5)
hold on
plot(date, Xtot(2,:), LineWidth=1.5)
ylabel('$\char"0023$ of cases','Interpreter','latex')
title('\textbf{Curve Fitting - Infected}','Interpreter','latex')
grid on
legend('Real Data', 'EKF Estimated Data','Interpreter','latex', 'Location','northeast')
xlim([date(1), date(end)])
set(gca, 'TickLabelInterpreter', 'Latex')


% Diagnosed - D
figure(3)
plot(date, D_data, LineWidth=1.5)
hold on
plot(date, Xtot(3,:), LineWidth=1.5)
ylabel('$\char"0023$ of cases','Interpreter','latex')
title('\textbf{Curve Fitting - Diagnosed}','Interpreter','latex')
grid on
legend('Real Data', 'EKF Estimated Data','Interpreter','latex', 'Location','northeast')
xlim([date(1), date(end)])
set(gca, 'TickLabelInterpreter', 'Latex')


% Lightly Threatned - T1
figure(4)
plot(date, T1_data, LineWidth=1.5)
hold on
plot(date, Xtot(4,:), LineWidth=1.5)
ylabel('$\char"0023$ of cases','Interpreter','latex')
title('\textbf{Curve Fitting - Hospitalised }','Interpreter','latex')
grid on
legend('Real Data', 'EKF Estimated Data','Interpreter','latex', 'Location','northeast')
xlim([date(1), date(end)])
set(gca, 'TickLabelInterpreter', 'Latex')


% ICU / Heavily Threatned - T2
figure(5)
plot(date, T2_data, LineWidth=1.5)
hold on
plot(date, Xtot(5,:), LineWidth=1.5)
ylabel('$\char"0023$ of cases','Interpreter','latex')
title('\textbf{Curve Fitting - ICUs / Heavily Threatned }','Interpreter','latex')
grid on
legend('Real Data', 'EKF Estimated Data','Interpreter','latex', 'Location','northeast')
xlim([date(1), date(end)])
set(gca, 'TickLabelInterpreter', 'Latex')


% Healed - H
figure(6)
plot(date, H_data, LineWidth=1.5)
hold on
plot(date, Xtot(6,:), LineWidth=1.5)
ylabel('$\char"0023$ of cases','Interpreter','latex')
title('\textbf{Curve Fitting - Healed}','Interpreter','latex')
grid on
legend('Real Data', 'EKF Estimated Data','Interpreter','latex', 'Location','northeast')
xlim([date(1), date(end)])
set(gca, 'TickLabelInterpreter', 'Latex')


% Deceased - E
figure(7)
plot(date, E_data, LineWidth=1.5)
hold on
plot(date, Xtot(7,:), LineWidth=1.5)
ylabel('$\char"0023$ of cases','Interpreter','latex')
title('\textbf{Curve Fitting - Expired}','Interpreter','latex')
grid on
legend('Real Data', 'EKF Estimated Data','Interpreter','latex', 'Location','northeast')
xlim([date(1), date(end)])
set(gca, 'TickLabelInterpreter', 'Latex')


% Coefficients alpha
figure(7)
plot(date, Xtot(8,:), LineWidth=1.5)
ylabel('Coefficients Values','Interpreter','latex')
title('\textbf{Coefficient $\alpha$}','Interpreter','latex')
grid on
legend('$\alpha$','Interpreter','latex', 'Location','northeast')
xlim([date(1), date(end)])
set(gca, 'TickLabelInterpreter', 'Latex')


% Coefficients beta
figure(8)
plot(date, Xtot(9,:), LineWidth=1.5)
ylabel('Coefficients Values','Interpreter','latex')
title('\textbf{Coefficient $\beta$}','Interpreter','latex')
grid on
legend('$\beta$','Interpreter','latex', 'Location','northeast')
xlim([date(1), date(end)])
set(gca, 'TickLabelInterpreter', 'Latex')


% Coefficient gamma
figure(9)
plot(date, Xtot(10,:), LineWidth=1.5)
ylabel('Coefficients Values','Interpreter','latex')
title('\textbf{Coefficient $\gamma$}','Interpreter','latex')
grid on
legend('$\gamma$','Interpreter','latex', 'Location','northeast')
xlim([date(1), date(end)])
set(gca, 'TickLabelInterpreter', 'Latex')


% Coefficients delta1 and delta2
figure(10)
plot(date, Xtot(11,:), LineWidth=1.5)
hold on
plot(date, Xtot(12,:), LineWidth=1.5)
ylabel('Coefficients Values','Interpreter','latex')
title('\textbf{Coefficients $\delta_1$ and $\delta_2$}','Interpreter','latex')
grid on
legend('$\delta_1$', '$\delta_2$','Interpreter','latex', 'Location','northeast')
xlim([date(1), date(end)])
set(gca, 'TickLabelInterpreter', 'Latex')


% Coefficients sigma1 and sigma2
figure(12)
plot(date, Xtot(14,:), LineWidth=1.5)
hold on
plot(date, Xtot(15,:), LineWidth=1.5)
ylabel('Coefficients Values','Interpreter','latex')
title('\textbf{Coefficients $\sigma_1$ and $\sigma_2$}','Interpreter','latex')
grid on
legend('$\sigma_1$','$\sigma_2$','Interpreter','latex', 'Location','northeast')
xlim([date(1), date(end)])
set(gca, 'TickLabelInterpreter', 'Latex')


% Coefficients tau1 and tau2
figure(13)
plot(date, Xtot(16,:), LineWidth=1.5)
hold on
plot(date, Xtot(17,:), LineWidth=1.5)
ylabel('Coefficients Values','Interpreter','latex')
title('\textbf{Coefficients $\tau_1$ and $\tau_2$}','Interpreter','latex')
grid on
legend('$\tau_1$','$\tau_2$','Interpreter','latex', 'Location','northeast')
xlim([date(1), date(end)])
set(gca, 'TickLabelInterpreter', 'Latex')

% Coefficients alpha
figure(14)
plot(date, Xtot(14,:), LineWidth=1.5)
ylabel('Coefficients Values','Interpreter','latex')
title('\textbf{Coefficient $\epsilon$}','Interpreter','latex')
grid on
legend('$\epsilon$','Interpreter','latex', 'Location','northeast')
xlim([date(1), date(end)])
set(gca, 'TickLabelInterpreter', 'Latex')
