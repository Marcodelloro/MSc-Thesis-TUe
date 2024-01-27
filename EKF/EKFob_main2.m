%% EKF - estimation of unknown parameters
% Improved code for EKF to estimate unknown parameters

clear all
close all
clc

load('/Users/marcodelloro/Desktop/Thesis/MSc-Thesis-TUe/Data_Collection/SIDTTHE_data.mat');
load('/Users/marcodelloro/Desktop/Thesis/MSc-Thesis-TUe/Data_Collection/Tests_data.mat');
load('/Users/marcodelloro/Desktop/Thesis/MSc-Thesis-TUe/Model_Fitting/EKF_coefs.mat');

addpath('/Users/marcodelloro/Downloads/casadi-3.6.3-osx64-matlab2018b')

date = SIDTTHE_data{1,1}.date;
Npop = 59240329; % Total Population of Italy

I_data = SIDTTHE_data{1,1}.data / Npop;      % I group - Infected Population
D_data = SIDTTHE_data{2,1}.data / Npop;      % D group - Diagnosed Population
T1_data = SIDTTHE_data{3,1}.data / Npop;     % T1 group - Hospitalised (Lightly Threatned) Population
T2_data = SIDTTHE_data{4,1}.data / Npop;     % T2 group - ICUs (Severly Threatned) Population
H_data = SIDTTHE_data{6,1}.data / Npop;      % T2 group - Healed Population
E_data = SIDTTHE_data{5,1}.data  / Npop;     % E group - Expired (Deceased) Population
S_data = ones(length(I_data),1)' - (I_data + D_data + T1_data + T2_data + H_data + E_data );   % S group - Susceptible Population

%% Normalisation of data

I_data = I_data./max(I_data);
D_data = D_data./max(D_data);
T1_data = T1_data./max(T1_data);
T2_data = T2_data./max(T2_data);
H_data = H_data./max(H_data);
E_data = E_data./max(E_data);
S_data = S_data./max(S_data);

%% Symbolic Equations 
% Define the set of differential equations and relative syms variables

syms S I D T1 T2 H E % Variables
syms alpha beta gamma lambda delta1 delta2 sigma1 tau1 epsilon sigma2 tau2 % Coefs

dS = -S * (alpha * I + beta * D);
dI = S * (alpha * I + beta * D) - (gamma + lambda) * I;
dD = I * gamma - D * (lambda + delta1 + delta2);
dT1 = delta1 * D - (sigma1 + tau1 + epsilon) * T1;
dT2 = delta2 * D - (sigma2 + tau2) * T2 + epsilon * T1;
dH = (I + D) * lambda + T1 * sigma1 + T2 * sigma2;
dE = T1 * tau1 + T2 * tau2;

% Differential equations for the parameters, which are zero
dAlpha = 0;
dGamma = 0;

% vetors for the computation - This Select what you want to compute

eqns = [dS; dI; dD; dT1; dT2; dH; dE; dAlpha; dGamma];
vars = [S I D T1 T2 H E alpha gamma];
static_coefs = [beta delta1 delta2 epsilon sigma1 sigma2 tau1 tau2 lambda]; % Coefs

%% Algorithm initialization

nstates = 9; % number of states - already augmented
nmeasures = 7; % number of measurements - already augmented
nparams = 2; % number of parameters that are estimated by the KF 

% States and Covariance initialization
x0A = [S_data(1); I_data(1); D_data(1); T1_data(1); T2_data(1); H_data(1); E_data(1); 0.6; 0.1];
P0 = diag([5 50 5 5 5 5 5 50 50]);

% Noise & Disturbances on the States
n = 0.0001; % Noise amplitude in the measurements - measurement without noise

w_alpha = 0.001; % Artificial disturbances on dynamics of coefficients
w_gamma = 0.01;

%% Computation of all Matrices & Vectors needed for EKF

w_dist = [w_alpha w_gamma]; % vector of the process disturbances

Zmeas = [S_data; I_data; D_data; T1_data; T2_data; H_data; E_data];

J = jacobian(eqns, vars);

L = [zeros((nstates - size(w_dist,2)),nparams); eye(nparams) ];

R = n*eye(nmeasures);  % Artificial disturbances on alpha dynamics

Qp = diag([w_dist]);   % Process noise covariance

C = [eye(nmeasures) zeros(nmeasures, nparams)];

%% Comptutaion of the EKF

% static parameters that are guessed 

beta_val = 0.5;
delta1_val = 0.002;
delta2_val = 5e-4;
epsi_val = 5e-4;
sigma1_val = 0.03;
sigma2_val = 0.01;
tau1_val = 0.03;
tau2_val = 0.008;
lambda_val = 0.05;

guessed_vals = [beta_val, delta1_val, delta2_val, epsi_val,...
                sigma1_val, sigma2_val, tau1_val, tau2_val, lambda_val];

for ii = 1:399

    Zmeas_k = Zmeas(:,ii);

    [Xk_plus, Pk_plus,Ob] = sidttheEKFobsv(x0A,P0,R,J,L,Qp,C,Zmeas_k,vars,static_coefs,guessed_vals,nstates);

    x0A = Xk_plus;
    P0 = Pk_plus;

    Xtot(:,ii) = Xk_plus;

end

%% Observability check

if rank(Ob) == nstates
    str = ['System is OBSERVABLE, rank = ', num2str(rank(Ob))];
    disp(str)
else
    str = ['ATTENTION !! System is NOT OBSERVABLE, rank = ', num2str(rank(Ob))];
    disp(str)
end

