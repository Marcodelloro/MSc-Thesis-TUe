%% DATA FROM ITALY
% The following dataset is referred to the COVID-19 pandemic occourred in
% Italy. The dataset has been retrived from https://lab24.ilsole24ore.com/coronavirus/?refresh_ce
% The data will be used to fit a simplification of the SIDARTHE model from the G.Giordano paper.

% Our interest are the 57 weeks between 31 august 2020 and 3 october 2021, so the year of
% COVID-19 outbrake, taking into account 2021 second wave and intial vaccination of the population.

clc
clear all 
close all

dataset = readtable('/Users/marcodelloro/Desktop/Thesis/Codes/Italy_complete_dataset.xlsx');
dataset = dataset(190:588,:);
dataset.data = datetime(dataset.data, 'InputFormat', 'dd-MMM-yyyy HH:mm:ss');

% Averaging of the data by week
window = 7; % weekly average

for ii = 1:window:size(dataset,1)
    hosp_avg.data(ii)= mean(dataset.ricoverati(ii:ii + window -1));             % weekly average of hospitalised pop
    ICU_avg.data(ii) = mean(dataset.terapia_intensiva(ii:ii+ window-1));        % weekly average of ICU pop
    pos_avg.data(ii) = mean(dataset.isolamento_domiciliare(ii:ii+ window-1));   % weekly average of Positive, detected, NOT THREATNED pop
    healed_avg.data(ii) = mean(dataset.guariti(ii:ii+ window-1));               % weekly average of healed pop
    dec_avg.data(ii) = mean(dataset.deceduti(ii:ii+ window-1));                 % weekly average of deceased pop
end

hosp_avg.data = nonzeros(hosp_avg.data);
ICU_avg.data = nonzeros(ICU_avg.data);
pos_avg.data = nonzeros(pos_avg.data);
healed_avg.data = nonzeros(healed_avg.data);
dec_avg.data = nonzeros(dec_avg.data);

t1 = datetime(2020,8,31);
t2 = datetime(2021,10,3);

hosp_avg.date = t1:caldays(7):t2;
ICU_avg.date = t1:caldays(7):t2;
pos_avg.date = t1:caldays(7):t2;
healed_avg.date = t1:caldays(7):t2;
dec_avg.date = t1:caldays(7):t2;


%% Estimate of the total number of infected, undetected (Not tested)

% This part of code is based on the paper "Estimating the undetected infections in the Covid-19 outbreak by harnessing capture–recapture methods"
% by Dankmar Böhning et Al.
% The Following code purpose is estimate data for "I" (Infected, undetected, asymptomatic, capable of infecting), in the Netherlands.
% For more info on the capture–recapture (CR) methods, refer to the original article cited above

% Initial data needed 
% Cumulative counts of infections = dec_avg.csum
% Cumulative counts of deceased = dec_avg.csum
% delta_N(t) = N(t) - N(t-1) -- New cases per day (t)
% delta_D(t) = D(t) - D(t-1) -- Deceased per day (t)

% H(t) =  [ delta_N(t) * ( delta_N(t) -1 ) ] / [ 1 + delta_N(t-1) - delta_D(t)] 
% H(t) = Hidden cases bias-corrected form by Chao

delta_N = dataset.nuovi_positivi;
delta_D = dataset.diff_deceduti;

H_start = ( delta_N(1) * (delta_N(1) - 1) ) / ( 1 + delta_N(1) - delta_D(1) );

varH_start = ( delta_N(1)^4 ) / ( ( 1 + delta_N(1) - delta_D(1) )^3)  + ...
               ( 4 * delta_N(1)^3 ) / ( ( 1 + delta_N(1) - delta_D(1) )^2) + ...
               ( delta_N(1)^2 ) / ( ( 1 + delta_N(1) - delta_D(1) ) );

% Note that also variance for H must be taken into account, so following
% the paper guidlines, variance estimate proposed by Niwitpong

% varH(t) = ( delta_N(t)^4 ) / ( 1 + delta_N(t-1) - delta_D(t) )^3  + ...
%           ( 4 * delta_N(t)^3 ) / ( 1 + delta_N(t-1) - delta_D(t) )^2 + ...
%           ( delta_N(t)^2 ) / ( 1 + delta_N(t-1) - delta_D(t) );

for ii = 2:size(delta_N,1)
    factor = delta_N(ii-1) - delta_D(ii);
        if factor < 0
           factor = 0;
        end
    % Hidden cases 

    num_H = delta_N(ii) * (delta_N(ii) - 1);
    den_H = 1 + factor;
    H(ii) = num_H/den_H;

    % Hidden cases variance

    varH(ii) = ( delta_N(ii)^4 ) / (den_H^3)  + ...
               ( 4 * delta_N(ii)^3 ) / (den_H^2) + ...
               ( delta_N(ii)^2 ) / (den_H);

    H_plus(ii) = H(ii) + 2.576*sqrt(varH(ii)); % 99perc confidence interval
    H_minus(ii) = H(ii) - 2.576*sqrt(varH(ii));

end

H(1) = H_start;
varH(1) = varH_start;
H_plus(1) = H(1) + 2.576*sqrt(varH(1)); % 99perc confidence interval
H_minus(1) = H(1) - 2.576*sqrt(varH(1)); 

% Averaging "Dark Number Data" - Hidden cases per week 

for ii = 1:window:size(dataset,1)
    dark_avg.data(ii) = mean(H(ii:ii + window -1));             % Weekly averaged Undetected Positive
    dark_plus_avg.data(ii) = mean(H_plus(ii:ii + window -1));   % Weekly averaged Undetected Positive - positive bound
    dark_minus_avg.data(ii) = mean(H_minus(ii:ii + window -1)); % Weekly averaged Undetected Positive - negative bound
    delta_N_avg.data(ii) = mean(delta_N(ii:ii + window -1));
end

dark_avg.data = nonzeros(dark_avg.data);                        % Number of HIDDEN Population PER DAY
dark_plus_avg.data = nonzeros(dark_plus_avg.data);
dark_minus_avg.data = nonzeros(dark_minus_avg.data);
delta_N_avg.data = nonzeros(delta_N_avg.data);

% Following the 'Linear Correlation approach' i calculate the ODE for infected
ratio_observed = (dark_avg.data + delta_N_avg.data) ./ delta_N_avg.data;
ratio_undetected = dark_avg.data ./ delta_N_avg.data;

infected_avg.date = pos_avg.date;
infected_avg.data = pos_avg.data .* ratio_undetected;

% We are going to assume a ratio also between people that are healed and
% diagnosed as the same rate for the infected non diaggnosed

heal_rate = healed_avg.data./pos_avg.data;
heal_notdiag = heal_rate .* infected_avg.data;

total_heal.data = healed_avg.data + heal_notdiag;
total_heal.date = healed_avg.date;

%% Saving all data

SIDTTHE_data = [infected_avg, pos_avg, hosp_avg, ICU_avg, dec_avg, total_heal];
filename = 'SIDTTHE_data.mat';
save(filename, 'SIDTTHE_data'); 

%% Plot section of the different trends

% "D" Positive, detected, NOT THREATNED population
figure()
plot(dataset.data, dataset.isolamento_domiciliare, LineWidth=1.5)
hold on 
plot(pos_avg.date, pos_avg.data, LineWidth=1.5)
ylabel('$\char"0023$ of cases', 'Interpreter', 'latex')
title('\textbf{Positive - Detected Population}', 'Interpreter', 'latex')
grid on
legend('Total', 'Averaged','Interpreter', 'latex')
xlim([pos_avg.date(1), pos_avg.date(end)])
set(gca, 'TickLabelInterpreter', 'Latex')


% "H" Healed population
figure()
plot(dataset.data, dataset.guariti, LineWidth=1.5)
hold on 
plot(healed_avg.date, healed_avg.data, LineWidth=1.5)
ylabel('$\char"0023$ of cases', 'Interpreter', 'latex')
title('\textbf{Healed Population}', 'Interpreter', 'latex')
grid on
legend('Total', 'Averaged','Interpreter', 'latex', 'Location','southeast')
xlim([pos_avg.date(1), pos_avg.date(end)])
set(gca, 'TickLabelInterpreter', 'Latex')


% "E" Expired (Deceased) population
figure()
plot(dataset.data, dataset.deceduti, LineWidth=1.5)
hold on 
plot(dec_avg.date, dec_avg.data, LineWidth=1.5)
ylabel('$\char"0023$ of cases', 'Interpreter', 'latex')
title('\textbf{Healed Population}', 'Interpreter', 'latex')
grid on
legend('Total', 'Averaged','Interpreter', 'latex', 'Location','southeast')
xlim([pos_avg.date(1), pos_avg.date(end)])
set(gca, 'TickLabelInterpreter', 'Latex')


% "T1" Hospitalised population
figure()
plot(dataset.data, dataset.ricoverati, LineWidth=1.5)
hold on 
plot(hosp_avg.date, hosp_avg.data, LineWidth=1.5)
ylabel('$\char"0023$ of cases', 'Interpreter', 'latex')
title('\textbf{Hospitalised Population}', 'Interpreter', 'latex')
grid on
legend('Total', 'Averaged','Interpreter', 'latex', 'Location','southeast')
xlim([pos_avg.date(1), pos_avg.date(end)])
set(gca, 'TickLabelInterpreter', 'Latex')


% "T2" Hospitalised population
figure()
plot(dataset.data, dataset.terapia_intensiva, LineWidth=1.5)
hold on 
plot(ICU_avg.date, ICU_avg.data, LineWidth=1.5)
ylabel('$\char"0023$ of cases', 'Interpreter', 'latex')
title('\textbf{Intensive Care Population}', 'Interpreter', 'latex')
grid on
legend('Total', 'Averaged','Interpreter', 'latex')
xlim([pos_avg.date(1), pos_avg.date(end)])
set(gca, 'TickLabelInterpreter', 'Latex')


% "I" Undetected infected, estimated model

var_area.x = [pos_avg.date flip(pos_avg.date)];
var_area.y = [dark_plus_avg.data' flip(dark_minus_avg.data')];

figure()
fill(var_area.x,var_area.y,[0 0.4470 0.7410],'FaceAlpha',.5,'EdgeColor', 'none','HandleVisibility', 'off')
hold on
plot(pos_avg.date, dark_avg.data, LineWidth=1.5, Color=[0 0.4470 0.7410])
ylabel('$\char"0023$ of cases','Interpreter','latex')
title('\textbf{Estimated Hidden Cases per day}','Interpreter','latex')
grid on
legend('Hidden Cases','Interpreter','latex','Location','northeast')
xlim([pos_avg.date(1), pos_avg.date(end)])
set(gca, 'TickLabelInterpreter', 'Latex')

% Ratio between "ESTIMATED TOTAL INFECTED / DETECTED INFECTED" - Trend

figure()
plot(pos_avg.date, ratio_observed, LineWidth=1.5, Color=[0 0.4470 0.7410])
ylabel('Ratio','Interpreter','latex')
title('\textbf{Ratio Estimated-Detected Infections}','Interpreter','latex')
grid on
legend('Ratio','Interpreter','latex')
xlim([pos_avg.date(1), pos_avg.date(end)])
set(gca, 'TickLabelInterpreter', 'Latex')


% Dark Numbers vs Actual Positive Trend

figure()
plot(pos_avg.date, delta_N_avg.data, LineWidth=1.5, Color=[0 0.4470 0.7410])
hold on
plot(pos_avg.date, dark_avg.data, LineWidth=1.5)
ylabel('Ratio','Interpreter','latex')
title('\textbf{Comparison Dark Data - Real Data}','Interpreter','latex')
grid on
legend('Real Data', 'Dark Data','Interpreter','latex')
xlim([pos_avg.date(1), pos_avg.date(end)])
set(gca, 'TickLabelInterpreter', 'Latex')

% Plot comparison between detected/diagnosed and infected population

figure()
plot(pos_avg.date, pos_avg.data, LineWidth=1.5, Color=[0 0.4470 0.7410])
hold on
plot(pos_avg.date, infected_avg.data, LineWidth=1.5)
ylabel('$\char"0023$ of cases','Interpreter','latex')
title('\textbf{Comparison Dark Data - Real Data}','Interpreter','latex')
grid on
legend('Real Data - Diagnosed', 'Dark Data - Undetected','Interpreter','latex')
xlim([pos_avg.date(1), pos_avg.date(end)])
set(gca, 'TickLabelInterpreter', 'Latex')

% Plot comparison between healed and "new healed" (based on assumed calculation)

figure()
plot(healed_avg.date, healed_avg.data, LineWidth=1.5, Color=[0 0.4470 0.7410])
hold on
plot(healed_avg.date, heal_notdiag, LineWidth=1.5)
ylabel('$\char"0023$ of cases','Interpreter','latex')
title('\textbf{Comparison Dark Data - Real Data}','Interpreter','latex')
grid on
legend('Real Data - Diagnosed', 'Dark Data - Undetected','Interpreter','latex', 'Location','southeast')
xlim([pos_avg.date(1), pos_avg.date(end)])
set(gca, 'TickLabelInterpreter', 'Latex')






