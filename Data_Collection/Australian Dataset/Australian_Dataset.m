clc 
clear all
close all

%% AUSTRALIAN DATA 

opts = detectImportOptions('/Users/marcodelloro/Desktop/National.csv');
opts = setvartype(opts, 'char'); % Set all variables to import as char
dataAUS = readtable('/Users/marcodelloro/Desktop/National.csv', opts);
dataAUS.Properties.VariableNames{'Var1'} = 'Date';
hospAUS = readtable('/Users/marcodelloro/Desktop/Hosp_Aus.csv', opts);
icuAUS =  readtable('/Users/marcodelloro/Desktop/ICU_Aus.csv', opts);

dataAUS.Hospitalised = hospAUS.AUS;
dataAUS.ICUs = icuAUS.AUS;

%Fixing the NaN values
columnNames = {'ActiveCases', 'Deaths', 'Recoveries','InactivesRemoved','Hospitalised','ICUs'}; % Replace with your actual column names
for i = 1:length(columnNames)
    columnName = columnNames{i};
    dataAUS.(columnName) = str2double(strrep(dataAUS.(columnName), ',', ''));
end

% Modifications of first "dataAUS" dataset
startDate = datetime(2020, 3, 15);
timeVector = startDate + days(0:size(dataAUS,1)-1); % Daily frequency
dataAUS.Date= timeVector';
dataAUS = dataAUS(dataAUS.Date >= datetime(2020, 8, 31) & dataAUS.Date <= datetime(2021, 10, 3), :);

dataAUS.Deaths = cumsum(dataAUS.Deaths);
dataAUS.Recoveries = cumsum(dataAUS.Recoveries);

DataStructAUS.Raw = dataAUS; %Saving Raw data  

%% Filtering - Smoothing of the data

% Cubic smoothing 
p = 0.4;
xi = [0:1:398];

dataAUS.ActiveCases = csaps(xi,log(dataAUS.ActiveCases),p,xi)';
dataAUS.Deaths = csaps(xi,log(dataAUS.Deaths),p,xi)';
dataAUS.Recoveries = csaps(xi,log(dataAUS.Recoveries),p,xi)';
dataAUS.Hospitalised = csaps(xi,log(dataAUS.Hospitalised),p,xi)';
dataAUS.ICUs = csaps(xi,log(dataAUS.ICUs),p,xi)';

% Apply a Savitzky - Golay filter on the cubic data
order = 5;
framelen = 27;

dataAUS.ActiveCases = sgolayfilt(dataAUS.ActiveCases,order,framelen);
dataAUS.Deaths = sgolayfilt(dataAUS.Deaths,order,framelen);
dataAUS.Recoveries = sgolayfilt(dataAUS.Recoveries,order,framelen);
dataAUS.Hospitalised = sgolayfilt(dataAUS.Hospitalised,order,framelen);
dataAUS.ICUs = sgolayfilt(dataAUS.ICUs,order,framelen);

dataAUS.ActiveCases = exp(dataAUS.ActiveCases);
dataAUS.Deaths = exp(dataAUS.Deaths);
dataAUS.Recoveries = exp(dataAUS.Recoveries);
dataAUS.Hospitalised = exp(dataAUS.Hospitalised);
dataAUS.ICUs = exp(dataAUS.ICUs);

DataStructAUS.Filterd = dataAUS; %Saving Filtered data  

%% Plots

% "D" Positive, detected, NOT THREATNED population
figure()
plot(dataAUS.Date, dataAUS.ActiveCases, LineWidth=1.5)
ylabel('$\char"0023$ of cases', 'Interpreter', 'latex')
title('\textbf{Positive - Detected Population}', 'Interpreter', 'latex')
grid on
xlim([dataAUS.Date(1), dataAUS.Date(end)])
set(gca, 'TickLabelInterpreter', 'Latex')


% "H" Healed population
figure()
plot(dataAUS.Date, dataAUS.Recoveries, LineWidth=1.5)
ylabel('$\char"0023$ of cases', 'Interpreter', 'latex')
title('\textbf{Healed Population}', 'Interpreter', 'latex')
grid on
xlim([dataAUS.Date(1), dataAUS.Date(end)])
set(gca, 'TickLabelInterpreter', 'Latex')


% "E" Expired (Deceased) population
figure()
plot(dataAUS.Date, dataAUS.Deaths, LineWidth=1.5)
ylabel('$\char"0023$ of cases', 'Interpreter', 'latex')
title('\textbf{Deceased Population}', 'Interpreter', 'latex')
grid on
xlim([dataAUS.Date(1), dataAUS.Date(end)])
set(gca, 'TickLabelInterpreter', 'Latex')


% "T1" Hospitalised population
figure()
plot(dataAUS.Date, dataAUS.Hospitalised, LineWidth=1.5)
ylabel('$\char"0023$ of cases', 'Interpreter', 'latex')
title('\textbf{Hospitalised Population}', 'Interpreter', 'latex')
grid on
xlim([dataAUS.Date(1), dataAUS.Date(end)])
set(gca, 'TickLabelInterpreter', 'Latex')


% "T2" ICU population
figure()
plot(dataAUS.Date, dataAUS.ICUs, LineWidth=1.5)
ylabel('$\char"0023$ of cases', 'Interpreter', 'latex')
title('\textbf{Intensive Care Population}', 'Interpreter', 'latex')
grid on
xlim([dataAUS.Date(1), dataAUS.Date(end)])
set(gca, 'TickLabelInterpreter', 'Latex')
