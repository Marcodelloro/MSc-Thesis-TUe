%% TRIAL TO MAKE THE MODEL SWITCH WHEN THE STATE RECHES A THRESHOLD

clc
clear all
close all

%% Previous data loading of the curves that must be fitted
import casadi.*

load('/Users/marcodelloro/Desktop/Thesis/MSc-Thesis-TUe/Data_Colection/SIDTTHE_data.mat');
load('/Users/marcodelloro/Desktop/Thesis/MSc-Thesis-TUe/Data_Colection/Tests_data.mat');

N = 399;
ICUs = 2000 ;

% Normalisation of loaded data so they add up to 1

date = SIDTTHE_data{1,1}.date;
Npop = 59240329; % Total Population of Italy
ICUs = 2000/Npop; % Total number of ICUs in Italy @ 9 Ottobre 2020

I_data = SIDTTHE_data{1,1}.data / Npop;     % I group - Infected Population **
D_data = SIDTTHE_data{2,1}.data / Npop;     % D group - Diagnosed Population
T1_data = SIDTTHE_data{3,1}.data / Npop;    % T1 group - Hospitalised (Lightly Threatned) Population
T2_data = SIDTTHE_data{4,1}.data / Npop;    % T2 group - ICUs (Severly Threatned) Population
H_data = SIDTTHE_data{6,1}.data / Npop;     % T2 group - Healed Population **
E_data = SIDTTHE_data{5,1}.data / Npop;     % E group - Expired (Deceased) Population

S_data = ones(length(I_data),1) - (I_data + D_data + T1_data + T2_data + H_data + E_data); 

% ** the groups I and H are based on strong assumptions and dark numbers estimates

DataArray = {S_data, I_data, D_data, T1_data, T2_data, H_data, E_data};

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
            tests = interp1(SIDTTHE_data{1, 1}.date, tests_avg, new_time);
        
        end

         data = vertcat(upsampled_DataArray{:});

    else
        disp('Invalid selection.');
    end
end

% Setting of the initial conditions for the solver

ODEs_0 = [  Npop;...       % S
            0;...           % I
            0;...           % D
            0;...           % T1
            0;...           % T2
            0;...           % H
            0;...           % E

         ];                % Initial values for ODEs

coefs_0 = [ 0.513;...       % alpha
            0.011;...       % beta
            0.2405;...      % gamma
            0.005;...       % delta1
            0.001;...       % delta2
            0.02;...        % epsi
            0.03;...        % sigma1
            0.065;...       % sigma2
            0.015;...       % tau1
            0.05;...        % tau2
            0.0596          % lambda

         ];                % Initial values for coefs

%% Settings

% Collocation settings
problem_options = NosnocProblemOptions();
solver_options = NosnocSolverOptions();
problem_options.n_s = 4;   % RK45
problem_options.irk_scheme = IRKSchemes.RADAU_IIA;

problem_options.use_fesd = 1;
problem_options.cross_comp_mode = 3;
solver_options.mpcc_mode = 'Scholtes_ineq';
solver_options.print_level = 2;
solver_options.s_elastic_max = 1e1;                    
solver_options.sigma_0 = 1;
solver_options.comp_tol = 1e-9;
solver_options.N_homotopy = 10;

%% Generate Model

model = model_switching(ODEs_0,coefs_0,ICUs);

%% Simulation settings

N_finite_elements = 1;  % Number of finite elements N_FE per control interval/stage
T_sim = N;             % Length of the control horizon T in the OCP
N_sim = N;             % This is the number of control intervals ns

problem_options.dcs_mode = 'Stewart';
% problem_options.dcs_mode = 'Step';

problem_options.T_sim = T_sim;
problem_options.N_finite_elements = N_finite_elements;
problem_options.N_sim = N_sim;
solver_options.use_previous_solution_as_initial_guess = 1;

model.f_q = (F - DataArray)^2;

%% Call FESD Integrator

integrator = NosnocIntegrator(model, problem_options, solver_options, [], []);
[results,stats] = integrator.solve();

%% Get variables into main workspace

unfold_struct(model,'base');
unfold_struct(problem_options,'base');
unfold_struct(results,'base');  
