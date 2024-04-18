%% DATA FROM ITALY
% The following dataset is referred to the COVID-19 pandemic occourred in
% Italy. The dataset has been retrived from https://lab24.ilsole24ore.com/coronavirus/?refresh_ce
% The data will be used to fit a simplification of the SIDARTHE model from the G.Giordano paper.

% Our interest are the 57 weeks between 31 august 2020 and 3 october 2021, so the year of
% COVID-19 outbrake, taking into account 2021 second wave and intial vaccination of the population.

clc
clear all 
close all
dataset = readtable('/Users/marcodelloro/Desktop/Thesis/MSc-Thesis-TUe/Data_Collection/Italian Dataset/Italy_complete_dataset.xlsx');
variants_data = readtable('/Users/marcodelloro/Desktop/Thesis/MSc-Thesis-TUe/Data_Collection/Italian Dataset/Italy_complete_dataset.xlsx', 'Sheet','Variants');
variants_data = variants_data(1:1368,:);

vax = load("vaxData.mat");

dataset = dataset(190:588,:);
dataset.data = datetime(dataset.data, 'InputFormat', 'dd-MMM-yyyy HH:mm:ss');

%% Smoothing and filtering 
% Cubic smoothing 
p = 0.4;
xi = [0:1:398];

hosp.data = csaps(xi,dataset.ricoverati,p,xi);
ICU.data = csaps(xi,dataset.terapia_intensiva,p,xi);
pos.data = csaps(xi,dataset.isolamento_domiciliare,p,xi);
healed.data = csaps(xi,dataset.guariti,p,xi);
dec.data = csaps(xi,dataset.deceduti,p,xi);
tests.data = csaps(xi,dataset.diff_tamponi,p,xi);

ICU_cube = ICU.data;

% Apply a Savitzky - Golay filter on the cubic data
order = 5;
framelen = 27;

hosp.data = sgolayfilt(hosp.data,order,framelen);
ICU.data = sgolayfilt(ICU.data,order,framelen);
pos.data = sgolayfilt(pos.data,order,framelen);
healed.data = sgolayfilt(healed.data,order,framelen);
dec.data = sgolayfilt(dec.data,order,framelen);
tests.data = sgolayfilt(tests.data,order,framelen);

t1 = datetime(2020,8,31);
t2 = datetime(2021,10,3);

hosp.date = t1:caldays(1):t2;
ICU.date = t1:caldays(1):t2;
pos.date = t1:caldays(1):t2;
healed.date = t1:caldays(1):t2;
dec.date = t1:caldays(1):t2;

indexOfInterest = 110:170;
figure()
plot(dataset.data, dataset.terapia_intensiva, LineWidth=2)
hold on 
plot(pos.date, ICU_cube, LineWidth=2)
hold on 
plot(pos.date, ICU.data, LineWidth=2)
ylabel('$\char"0023$ of cases', 'Interpreter', 'latex')
% title('\textbf{Intensive Care Units}', 'Interpreter', 'latex')
grid on
legend('Raw', 'Cubic Spline','S-G Filter','Interpreter', 'latex')
xlim([pos.date(1), pos.date(end)])
set(gca, 'TickLabelInterpreter', 'Latex')
insetAx = axes('position', [.35 .175 .25 .25]);
plot(insetAx, pos.date(indexOfInterest), dataset.terapia_intensiva(indexOfInterest), 'LineWidth', 2)
hold on
plot(insetAx, pos.date(indexOfInterest), ICU_cube(indexOfInterest), 'LineWidth', 2)
hold on
plot(insetAx, pos.date(indexOfInterest), ICU.data(indexOfInterest), 'LineWidth', 2)
axis tight
set(insetAx, 'XTick', [], 'YTick', []);  % Remove ticks

%% Estimate of the total number of infected, undetected (Not tested)

% This part of code is based on the paper "Estimating the undetected infections in the Covid-19 outbreak by harnessing capture–recapture methods"
% by Dankmar Böhning et Al.
% The Following code purpose is estimate data for "I" (Infected, undetected, asymptomatic, capable of infecting), in the Netherlands.
% For more info on the capture–recapture (CR) methods, refer to the original article cited above

% Initial data needed 
% Cumulative counts of infections = dec.csum
% Cumulative counts of deceased = dec.csum
% delta_N(t) = N(t) - N(t-1) -- New cases per day (t)
% delta_D(t) = D(t) - D(t-1) -- Deceased per day (t)

% H(t) =  [ delta_N(t) * ( delta_N(t) -1 ) ] / [ 1 + delta_N(t-1) - delta_D(t)] 
% H(t) = Hidden cases bias-corrected form by Chao

delta_N = dataset.nuovi_positivi;
delta_D = dataset.diff_deceduti;
delta_H = dataset.diff_guariti;

H_start = ( delta_N(1) * (delta_N(1) - 1) ) / ( 1 + delta_N(1) - delta_D(1) );

varH_start = ( delta_N(1)^4 ) / ( ( 1 + delta_N(1) - delta_D(1)  )^3)  + ...
               ( 4 * delta_N(1)^3 ) / ( ( 1 + delta_N(1) - delta_D(1) )^2) + ...
               ( delta_N(1)^2 ) / ( ( 1 + delta_N(1) - delta_D(1)  ) );

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

    num_H = delta_N(ii) * (delta_N(ii-1));
    den_H = 1 + factor;
    H(ii) = num_H/den_H;

    % Hidden cases variance

    varH(ii) = ( delta_N(ii)^4 ) / (den_H^3)  + ...
               ( 4 * delta_N(ii)^3 ) / (den_H^2) + ...
               ( delta_N(ii)^2 ) / (den_H);

    H_plus(ii) = H(ii) + 3.96*sqrt(varH(ii)); % 99perc confidence interval
    H_minus(ii) = H(ii) - 3.96*sqrt(varH(ii));

end

H(1) = H_start;
varH(1) = varH_start;
H_plus(1) = H(1) + 1.96*sqrt(varH(1)); % 99perc confidence interval
H_minus(1) = H(1) - 1.96*sqrt(varH(1)); 

figure(34)
scatter(dataset.data, H,40 ,"filled")
ylabel('$\textit{$\Delta$I}$', 'Interpreter', 'latex');
% title('$\textbf{New Undetected cases per day}$', 'Interpreter', 'latex');
xlim([dataset.data(1) dataset.data(end)])
ylim([min(H) 10e4])
set(gca, 'TickLabelInterpreter', 'Latex')
yLimits = ylim;
yticks(0:2e4:10e4);
grid on
box on

% Filtering "Dark Number Data" cubic spline + SGF - Hidden cases per week 
order2 = 3;
dark_avg.data = csaps(xi,H,p,xi);
dark_plus_avg.data = csaps(xi,H_plus,p,xi);
dark_minus_avg.data = csaps(xi,H_minus,p,xi);
delta_N_avg.data = csaps(xi,delta_N,p,xi);

dark_avg.data = sgolayfilt(dark_avg.data,order2,framelen);
dark_plus_avg.data = sgolayfilt(dark_plus_avg.data,order2,framelen);
dark_minus_avg.data = sgolayfilt(dark_minus_avg.data,order2,framelen);
delta_N_avg.data = sgolayfilt(delta_N_avg.data,order2,framelen);


% Following the 'Linear Correlation approach' i calculate the ODE for infected
ratio_observed = (dark_avg.data + delta_N_avg.data) ./ delta_N_avg.data;
ratio_undetected = dark_avg.data ./ delta_N_avg.data;

ratioPlus = dark_plus_avg.data./ delta_N_avg.data;
ratioMinus = dark_minus_avg.data./ delta_N_avg.data;

infected_avg.date = pos.date;
infected_avg.data = pos.data .* ratio_undetected;
infected_avg_plus = pos.data .* ratioPlus;
infected_avg_minus = pos.data .* ratioMinus;

figure()
scatter(dataset.data, dark_avg.data, 'filled')
ylabel('$\textit{$\Delta$H}$', 'Interpreter', 'latex');
title('$\textbf{Filtered \textit{New Undetected cases per day}$', 'Interpreter', 'latex');
xlim([dataset.data(1) dataset.data(end)])
set(gca, 'TickLabelInterpreter', 'Latex')
grid on

% We are going to assume a ratio also between people that are healed and
% diagnosed as the same rate for the infected non diagnosed

heal_rate = delta_H'./pos.data;
heal_notdiag = heal_rate .* infected_avg.data;
delta_H(1) = delta_H(1) + healed.data(1);
added_heal = round(delta_H + heal_notdiag');

total_heal.data = cumsum(added_heal);
total_heal.date = healed.date;

% Filtering the data so can be smoothed out 

total_heal.data = csaps(xi,total_heal.data,p,xi);
total_heal.data = sgolayfilt(total_heal.data,order,framelen);

figure()
plot(total_heal.date, total_heal.data, LineWidth=1.5)
hold on 
plot(pos.date, healed.data, LineWidth=1.5)
ylabel('$\char"0023$ of cases', 'Interpreter', 'latex')
title('\textbf{Observed vs Total Healed Population}', 'Interpreter', 'latex')
grid on
legend('Total', 'Observed','Interpreter', 'latex','location','southeast')
xlim([pos.date(1), pos.date(end)])
set(gca, 'TickLabelInterpreter', 'Latex')


%% Modifing and working on the variants

varNames = horzcat('date', variants_data.variant(1:23)');
variants = array2table(zeros(57, 24), 'VariableNames', varNames);
variants_data.percent_variant = str2double(variants_data.percent_variant);

variants.date = (datetime(2020,8,31):caldays(7):datetime(2021,10,3))';

idx = 1;

for ii=1:23:size(variants_data,1)-23

    perc = variants_data.percent_variant(ii:ii+22)';
    variants{idx,2:end} = perc;
    idx = idx+1;
end

% cancelling columns with no values

zeroColumns = [3 4 10 11 16 17 18 19 22 23 24];
variants(:, zeroColumns) = [];
variants =  variants(1:57,:);

%% Saving all data

SIDTTHE_data = {infected_avg, 'Infected'; pos, 'Diagnosed'; hosp, 'Hospitalised'; ICU, 'ICU threatend'; dec, 'Deceased'; total_heal, 'HealedAug'; healed, 'Healed'};
filename = 'SIDTTHE_data_DEF.mat';
save(filename, 'SIDTTHE_data'); 

% save of the data from Testing activities
filename = 'Tests_data.mat';
save(filename, 'tests'); 

%% Plot section of the different trends

% "D" Positive, detected, NOT THREATNED population
figure()
plot(dataset.data, dataset.isolamento_domiciliare, LineWidth=1.5)
hold on 
plot(pos.date, pos.data, LineWidth=1.5)
ylabel('$\char"0023$ of cases', 'Interpreter', 'latex')
title('\textbf{Positive - Detected Population}', 'Interpreter', 'latex')
grid on
legend('Total', 'Averaged','Interpreter', 'latex')
xlim([pos.date(1), pos.date(end)])
set(gca, 'TickLabelInterpreter', 'Latex')


% "H" Healed population
figure()
plot(dataset.data, dataset.guariti, LineWidth=1.5)
hold on 
plot(healed.date, healed.data, LineWidth=1.5)
ylabel('$\char"0023$ of cases', 'Interpreter', 'latex')
title('\textbf{Healed Population}', 'Interpreter', 'latex')
grid on
legend('Total', 'Averaged','Interpreter', 'latex', 'Location','southeast')
xlim([pos.date(1), pos.date(end)])
set(gca, 'TickLabelInterpreter', 'Latex')


% "E" Expired (Deceased) population
figure()
plot(dataset.data, dataset.deceduti, LineWidth=1.5)
hold on 
plot(dec.date, dec.data, LineWidth=1.5)
ylabel('$\char"0023$ of cases', 'Interpreter', 'latex')
title('\textbf{Healed Population}', 'Interpreter', 'latex')
grid on
legend('Total', 'Averaged','Interpreter', 'latex', 'Location','southeast')
xlim([pos.date(1), pos.date(end)])
set(gca, 'TickLabelInterpreter', 'Latex')


% "T1" Hospitalised population
figure()
plot(dataset.data, dataset.ricoverati, LineWidth=1.5)
hold on 
plot(hosp.date, hosp.data, LineWidth=1.5)
ylabel('$\char"0023$ of cases', 'Interpreter', 'latex')
title('\textbf{Hospitalised Population}', 'Interpreter', 'latex')
grid on
legend('Total', 'Averaged','Interpreter', 'latex', 'Location','southeast')
xlim([pos.date(1), pos.date(end)])
set(gca, 'TickLabelInterpreter', 'Latex')


% "T2" Hospitalised population
figure()
plot(dataset.data, dataset.terapia_intensiva, LineWidth=1.5)
hold on 
plot(ICU.date, ICU.data, LineWidth=1.5)
ylabel('$\char"0023$ of cases', 'Interpreter', 'latex')
title('\textbf{Intensive Care Population}', 'Interpreter', 'latex')
grid on
legend('Total', 'Averaged','Interpreter', 'latex')
xlim([pos.date(1), pos.date(end)])
set(gca, 'TickLabelInterpreter', 'Latex')


% "I" Undetected infected, estimated model

var_area.x = [pos.date flip(pos.date)];
var_area.y = [dark_plus_avg.data'; flip(dark_minus_avg.data')];
var_area.Ivar = [infected_avg_plus'; flip(infected_avg_minus')];

figure()
plot(pos.date, dark_avg.data, LineWidth=2, Color=[0 0.4470 0.7410])
hold on
fill(var_area.x,var_area.y,[0 0.4470 0.7410],'FaceAlpha',.3,'EdgeColor','none')
ylabel('$\char"0023$ of cases','Interpreter','latex')
% title('\textbf{Estimated Hidden Cases per day}','Interpreter','latex')
grid on
legend('Hidden Cases','95 \% Confidence Interval','Interpreter','latex','Location','northeast')
xlim([pos.date(1), pos.date(end)])
set(gca, 'TickLabelInterpreter', 'Latex')
box on

save('Var_Infected.mat', 'var_area');

% Ratio between "ESTIMATED TOTAL INFECTED / DETECTED INFECTED" - Trend

figure()
plot(pos.date, ratio_undetected, LineWidth=2, Color=[0 0.4470 0.7410])
% title('\textbf{Ratio $\zeta$ Undetected - Detected Infections }','Interpreter','latex')
grid on
yax = ylabel('{Ratio Level}','Interpreter','latex');
% legend('$\zeta$','Interpreter','latex')
xlim([pos.date(1), pos.date(end)])
set(gca, 'TickLabelInterpreter', 'Latex')


% Dark Numbers vs Actual Positive Trend

figure()
plot(pos.date, delta_N_avg.data, LineWidth=1.5, Color=[0 0.4470 0.7410])
hold on
plot(pos.date, dark_avg.data, LineWidth=1.5)
ylabel('Ratio','Interpreter','latex')
title('\textbf{Comparison Dark Data - Real Data}','Interpreter','latex')
grid on
legend('Real Data', 'Dark Data','Interpreter','latex')
xlim([pos.date(1), pos.date(end)])
set(gca, 'TickLabelInterpreter', 'Latex')

% Plot comparison between detected/diagnosed and infected population

figure(30)
plot(pos.date, pos.data,LineWidth=1.5, LineStyle="--",Color=[0 0 0])
hold on
plot(pos.date, infected_avg.data, Color=[0 0.4470 0.7410], LineWidth=1.5)
colorsct=[0 0.4470 0.7410];
hold on 
fill(var_area.x,var_area.Ivar, colorsct, 'FaceAlpha', .2, 'EdgeColor', 'none');
ylabel('$\char"0023$ of cases','Interpreter','latex')
grid on
legend('Real Data - Diagnosed', 'Dark Data - Undetected','95\% Confidence Interval','Interpreter','latex')
set(gca, 'TickLabelInterpreter', 'Latex')
xlim([pos.date(1), pos.date(end)])
box on

% Plot comparison between healed and "new healed" (based on assumed calculation)

figure()
plot(healed.date, healed.data, LineWidth=2, Color=[0 0.4470 0.7410])
hold on
plot(healed.date, total_heal.data, LineWidth=2)
ylabel('$\char"0023$ of cases','Interpreter','latex')
% title('\textbf{Observed vs Total Healed population}','Interpreter','latex')
grid on
legend('Observed Healed', 'Total Healed','Interpreter','latex', 'Location','southeast')
xlim([pos.date(1), pos.date(end)])
set(gca, 'TickLabelInterpreter', 'Latex')

%% Subplot for thesis

figure("Units","centimeters","Position",[0,0,30,40])
subplot(3,2,1)
plot(dataset.data, dataset.isolamento_domiciliare, LineWidth=2)
ylabel('$\char"0023$ of cases', 'Interpreter', 'latex')
title('\textbf{Positive - Detected Population}', 'Interpreter', 'latex')
grid on
xlim([pos.date(1), pos.date(end)])
set(gca, 'TickLabelInterpreter', 'Latex')

subplot(3,2,2)
plot(dataset.data, dataset.ricoverati, LineWidth=2)
ylabel('$\char"0023$ of cases', 'Interpreter', 'latex')
title('\textbf{Hospitalised Population}', 'Interpreter', 'latex')
grid on
xlim([pos.date(1), pos.date(end)])
set(gca, 'TickLabelInterpreter', 'Latex')

subplot(3,2,3)
plot(dataset.data, dataset.terapia_intensiva, LineWidth=2)
ylabel('$\char"0023$ of cases', 'Interpreter', 'latex')
title('\textbf{Intensive Care Population}', 'Interpreter', 'latex')
grid on
xlim([pos.date(1), pos.date(end)])
set(gca, 'TickLabelInterpreter', 'Latex')

subplot(3,2,4);
plot(dataset.data, dataset.guariti, LineWidth=2)
ylabel('$\char"0023$ of cases', 'Interpreter', 'latex')
title('\textbf{Healed Population}', 'Interpreter', 'latex')
grid on
xlim([pos.date(1), pos.date(end)])
set(gca, 'TickLabelInterpreter', 'Latex')

subplot(3,2,5);
plot(dataset.data, dataset.deceduti, LineWidth=2)
ylabel('$\char"0023$ of cases', 'Interpreter', 'latex')
title('\textbf{Deceased Population}', 'Interpreter', 'latex')
grid on
xlim([pos.date(1), pos.date(end)])
set(gca, 'TickLabelInterpreter', 'Latex')

subplot(3,2,6);
plot(dec.date, vax.vaxData.sum_d2, LineWidth=2)
hold on 
plot(dec.date, vax.vaxData.sum_dpi, LineWidth=1.5)
ylabel('$\char"0023$ of cases','Interpreter','latex')
title('\textbf{Vaccinated Population}','Interpreter','latex')
grid on
legend('$S$ Vaccines', '$H$ Vaccines','Interpreter','latex', 'Location','northwest')
xlim([pos.date(1), pos.date(end)])
set(gca, 'TickLabelInterpreter', 'Latex')
% sgtitle('\textbf{Raw Available Data}','Interpreter','latex');


%% Comparison plots between coefficients and other interesting parameters

% load("Optimization_Results.mat")
% Npop = 59240329; % Total Population of Italy
% 
% figure()
% plot(dataset.data, dataset.diff_tamponi, LineWidth=1.5)
% hold on
% plot(dataset.data, Opti_results{2,1}.gamma, LineWidth=1.5)
% hold on
% plot(healed.date, heal_notdiag, LineWidth=1.5)
% ylabel('$\char"0023$ of cases','Interpreter','latex')
% title('\textbf{Comparison Dark Data - Real Data}','Interpreter','latex')
% grid on
% legend('Real Data - Diagnosed', 'Dark Data - Undetected','Interpreter','latex', 'Location','southeast')
% xlim([pos.date(1), pos.date(end)])
% set(gca, 'TickLabelInterpreter', 'Latex')

%% SARS-CoV-V2 Variants figure


customColors = {[0,	0.447,	0.741],
                [0.85,	0.325,	0.098],
                [0.929,	0.694,	0.125],
                [0.494,	0.184,	0.556],
                [0.466,	0.674,	0.188],
                [0.301,	0.745,	0.933],
                [0.635,	0.078,	0.184],
                [0 0 1],            
                [1 0 0],              
                [1 1 0], 
                [0.1 0.1 0.1],
                [0 1 0.5]
                };


% modifiy data for the plot
data_barplot = removevars(variants,'date');
data_barplot = table2array(data_barplot);
date_barplot = variants.date;

figure()
h = bar(date_barplot, data_barplot, 0.8, "stacked");
for ii = 1:size(customColors, 1)
    set(h(1, ii), 'facecolor', customColors{ii, 1})
end
ylabel('$\%$ Variant', 'Interpreter', 'latex')
legend('SARS-CoV-2', 'B.1.351', 'B.1.617.2 - DELTA VARIANT', ...
    'B.1.1.7 - ALPHA VARIANT', 'B.1.427/B.1.429', 'B.1.525', ...
    'B.1.620', 'B.1.621', 'BA.1 - OMICRON BA1 VARIANT', ...
    'BA.2 - OMICRON BA2 VARIANT', 'C.37', 'P.1 - GAMMA VARIANT', ...
    'Interpreter', 'latex', 'Location', 'southoutside', 'NumColumns', 3)
xlim([variants.date(1), variants.date(end)])
ylim([0, 100.1])
set(gca, 'TickLabelInterpreter', 'Latex')

% Smoothing of Most importat variants to have a "Variants Plot"

% Cubic smoothing 

xii = [0:1:56];

SmoothVar.date = variants.date;
SmoothVar.SC2 = csaps(xii,variants.Other,p,xii);
SmoothVar.alpha = csaps(xii,variants.("B.1.1.7"),p,xii);
SmoothVar.gamma = csaps(xii,variants.("P.1"),p,xii);
SmoothVar.delta = csaps(xii,variants.("B.1.617.2"),p,xii);


save('SmoothVariantasIta.mat','SmoothVar');

%% avg OBS/TOT ratio

csumD_N = sum(delta_N);
csumD_H = sum(dark_avg.data);

r = (csumD_N + csumD_H) / csumD_N;