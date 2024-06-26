clc 
clear all
close all

set(0,'DefaultFigureWindowStyle','docked');

%% AUSTRALIAN DATA 

opts = detectImportOptions('/Users/marcodelloro/Desktop/Thesis/MSc-Thesis-TUe/Data_Collection/Australian Dataset/National.csv');
opts = setvartype(opts, 'char'); % Set all variables to import as char
dataAUS = readtable('/Users/marcodelloro/Desktop/Thesis/MSc-Thesis-TUe/Data_Collection/Australian Dataset/National.csv', opts);
dataAUS.Properties.VariableNames{'Var1'} = 'Date';
hospAUS = readtable('/Users/marcodelloro/Desktop/Thesis/MSc-Thesis-TUe/Data_Collection/Australian Dataset/Hosp_Aus.csv', opts);
icuAUS =  readtable('/Users/marcodelloro/Desktop/Thesis/MSc-Thesis-TUe/Data_Collection/Australian Dataset/ICU_Aus.csv', opts);
newpos = readtable('/Users/marcodelloro/Desktop/Thesis/MSc-Thesis-TUe/Data_Collection/Australian Dataset/PosNet.csv');
dataAUS.Hospitalised = hospAUS.AUS;
dataAUS.ICUs = icuAUS.AUS;

values = newpos.NewCases_Day;
dates2 = datetime('25-Jan-2020', 'InputFormat', 'dd-MMM-yyyy'):days(1):datetime('04-Aug-2023', 'InputFormat', 'dd-MMM-yyyy');
newpos = table(dates2', values, 'VariableNames', {'Date', 'NewPos'});
val2 = newpos.NewPos(newpos.Date >= datetime(2021, 07, 12) & newpos.Date <= datetime(2022, 08, 14), :);

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
dataAUS = dataAUS(dataAUS.Date >= datetime(2021, 07, 12) & dataAUS.Date <= datetime(2022, 08, 14), :);
newpos = table(dataAUS.Date, val2, 'VariableNames', {'Date', 'NewPos'});

% Fixing Healing data
Healed1 = dataAUS.InactivesRemoved(dataAUS.Date >= datetime(2021, 07, 12) & dataAUS.Date <= datetime(2022, 01, 17), :)+...
          dataAUS.Recoveries(dataAUS.Date >= datetime(2021, 07, 12) & dataAUS.Date <= datetime(2022, 01, 17), :);
Healed2 = dataAUS.InactivesRemoved(dataAUS.Date >= datetime(2022, 01, 18) & dataAUS.Date <= datetime(2022, 08, 14), :);

dataAUS.newHealed = [Healed1; Healed2];

DataStructAUS.Raw = dataAUS; %Saving Raw data  
%% Filtering - Smoothing of the data

% Cubic smoothing 
p = 0.4;
xi = [0:1:398];

dataAUS.ActiveCases = csaps(xi,log(dataAUS.ActiveCases),p,xi)';
dataAUS.Deaths = csaps(xi,log(dataAUS.Deaths),p,xi)';
dataAUS.Recoveries = csaps(xi,log(dataAUS.newHealed),p,xi)';
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

%% Vax data - Australian dataset

Vax_Data = readtable('/Users/marcodelloro/Desktop/Thesis/MSc-Thesis-TUe/Data_Collection/Australian Dataset/Vax2.csv', opts);
Vax_Data = Vax_Data(:, {'date', 'people_fully_vaccinated'});

% Modifications of "date" column
startDate2 = datetime(2021, 2, 21);
timeVector2 = startDate2 + days(0:size(Vax_Data,1)-1); % Daily frequency
Vax_Data.date= timeVector2';
Vax_Data = Vax_Data(Vax_Data.date >= datetime(2021, 07, 12) & Vax_Data.date <= datetime(2022, 08, 14), :);

Vax_Data.people_fully_vaccinated = str2double(strrep(Vax_Data.people_fully_vaccinated, ',', ''));

% Substitute NaN values
for i = 2:height(Vax_Data)  
    if isnan(Vax_Data{i, 'people_fully_vaccinated'})  % If the current value is NaN
        Vax_Data{i, 'people_fully_vaccinated'} = Vax_Data{i-1, 'people_fully_vaccinated'}; 
    end
end

% % Modifications of "date" column
% timeVector3 = datetime(2020, 8, 31):(Vax_Data.date(1)-1); % Daily frequency
% previousvax = zeros(length(timeVector3),1);
% prev_tab = table(timeVector3', previousvax, 'VariableNames', {'date', 'people_fully_vaccinated'});
% 
% % Full Vax Table
% Vax_Data = [prev_tab; Vax_Data];

%% Plots

% Double - Dosed Vaccines in Australia https://covidbaseau.com/
figure(45)
plot(Vax_Data.date, Vax_Data.people_fully_vaccinated, LineWidth=2)
ylabel('$\char"0023$ of cases', 'Interpreter', 'latex')
title('\textbf{Vaccinated Population}', 'Interpreter', 'latex')
grid on
xlim([Vax_Data.date(1), Vax_Data.date(end)])
set(gca, 'TickLabelInterpreter', 'Latex')

% Subplot for thesis

figure("Units","centimeters","Position",[0,0,30,40])
subplot(3,2,1)
plot(dataAUS.Date, dataAUS.ActiveCases, LineWidth=2)
ylabel('$\char"0023$ of cases', 'Interpreter', 'latex')
title('\textbf{Positive - Detected Population}', 'Interpreter', 'latex')
grid on
xlim([dataAUS.Date(1), dataAUS.Date(end)])
set(gca, 'TickLabelInterpreter', 'Latex')

subplot(3,2,2)
plot(dataAUS.Date, dataAUS.Hospitalised, LineWidth=2)
ylabel('$\char"0023$ of cases', 'Interpreter', 'latex')
title('\textbf{Hospitalised Population}', 'Interpreter', 'latex')
grid on
xlim([dataAUS.Date(1), dataAUS.Date(end)])
set(gca, 'TickLabelInterpreter', 'Latex')

subplot(3,2,3)
plot(dataAUS.Date, dataAUS.ICUs, LineWidth=2)
ylabel('$\char"0023$ of cases', 'Interpreter', 'latex')
title('\textbf{Intensive Care Population}', 'Interpreter', 'latex')
grid on
xlim([dataAUS.Date(1), dataAUS.Date(end)])
set(gca, 'TickLabelInterpreter', 'Latex')

subplot(3,2,4);
plot(dataAUS.Date, dataAUS.Recoveries, LineWidth=2)
ylabel('$\char"0023$ of cases', 'Interpreter', 'latex')
title('\textbf{Healed Population}', 'Interpreter', 'latex')
grid on
xlim([dataAUS.Date(1), dataAUS.Date(end)])
set(gca, 'TickLabelInterpreter', 'Latex')

subplot(3,2,5);
plot(dataAUS.Date, dataAUS.Deaths, LineWidth=2)
ylabel('$\char"0023$ of cases', 'Interpreter', 'latex')
title('\textbf{Deceased Population}', 'Interpreter', 'latex')
grid on
xlim([dataAUS.Date(1), dataAUS.Date(end)])
set(gca, 'TickLabelInterpreter', 'Latex')

subplot(3,2,6);
plot(Vax_Data.date, Vax_Data.people_fully_vaccinated, LineWidth=2)
ylabel('$\char"0023$ of cases','Interpreter','latex')
title('\textbf{Vaccinated Population}','Interpreter','latex')
grid on
xlim([dataAUS.Date(1), dataAUS.Date(end)])
set(gca, 'TickLabelInterpreter', 'Latex')
% sgtitle('\textbf{Raw Available Data}','Interpreter','latex');


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

delta_N = newpos.NewPos;
delta_D = diff(DataStructAUS.Raw.Deaths);
delta_D = [delta_D; 15];
delta_H = diff(DataStructAUS.Raw.newHealed);
delta_H = [delta_H; 879];

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

    H_plus(ii) = H(ii) + 1.96*sqrt(varH(ii)); % 99perc confidence interval
    H_minus(ii) = H(ii) - 1.96*sqrt(varH(ii));

end


figure(34)
scatter(DataStructAUS.Raw.Date, H,40 ,"filled",'MarkerFaceColor',[0.8500 0.3250 0.0980])
ylabel('$\textit{$\Delta$I}$', 'Interpreter', 'latex');
% title('$\textbf{New Undetected cases per day}$', 'Interpreter', 'latex');
xlim([DataStructAUS.Raw.Date(1) DataStructAUS.Raw.Date(end)])
ylim([0 1e5])
set(gca, 'TickLabelInterpreter', 'Latex')
grid on
box on
%% fixing numerical errors in Hplus/Hminus
H_plus(69) = H_plus(70);
H_plus(82) = H_plus(83);
H_plus(165) = H_plus(166);
 
H(1) = H_start;
varH(1) = varH_start;
H_plus(1) = H(1) + 1.96*sqrt(varH(1)); % 99perc confidence interval
H_minus(1) = H(1) - 1.96*sqrt(varH(1)); 

H_minus(H_minus < 0) = 0;

%% Plots
figure(10)
plot(dataAUS.Date, H, LineWidth=1.5)
ylabel('$\textit{$\Delta$I}$', 'Interpreter', 'latex');
title('$\textbf{New Undetected cases per day}$', 'Interpreter', 'latex');
xlim([dataAUS.Date(1), dataAUS.Date(end)])
grid on

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


% Following the 'Linear Correlation approach'
ratio_observed = (dark_avg.data + delta_N_avg.data) ./ delta_N_avg.data;
ratio_undetected = dark_avg.data ./ delta_N_avg.data;
ratio_undetected2 = dark_minus_avg.data ./ delta_N_avg.data;
ratioPlus = dark_plus_avg.data./ delta_N_avg.data;
ratioMinus = dark_minus_avg.data./ delta_N_avg.data;

infected_avg = DataStructAUS.Filterd.ActiveCases .* ratio_undetected';
infected_avg_plus = DataStructAUS.Filterd.ActiveCases .* ratioPlus'*1.05;
infected_avg_minus = DataStructAUS.Filterd.ActiveCases .* ratioMinus'*0.95;
% 
% infected_avg_plus(212:246) = infected_avg(212:246)*1.2;
% infected_avg_minus(1:160) = infected_avg(1:160)*0.5;
figure(11)
scatter(dataAUS.Date, dark_avg.data, 'filled')
ylabel('$\textit{$\Delta$H}$', 'Interpreter', 'latex');
yax = ylabel('\textbf{$\mathbf{\char"0023}$ of cases}','Interpreter','latex');
yax.FontSize = 14;
% title('\textbf{Filtered \textit{New Undetected cases per day}}', 'Interpreter', 'latex');
xlim([dataAUS.Date(1), dataAUS.Date(end)])
set(gca, 'TickLabelInterpreter', 'Latex')
grid on

% We are going to assume a ratio also between people that are healed and
% diagnosed as the same rate for the infected non diagnosed

heal_rate = delta_H./DataStructAUS.Filterd.ActiveCases;
heal_notdiag = heal_rate .* infected_avg;
delta_H(1) = delta_H(1) + DataStructAUS.Filterd.newHealed(1);
added_heal = round(delta_H + heal_notdiag);

total_heal.data = cumsum(added_heal);
total_heal.date = DataStructAUS.Filterd.Date;

% Filtering the data so can be smoothed out 

total_heal.data = csaps(xi,total_heal.data,p,xi);
total_heal.data = sgolayfilt(total_heal.data,order,framelen);

figure(12)
plot(dataAUS.Date, DataStructAUS.Filterd.Recoveries, LineWidth=2)
hold on
plot(dataAUS.Date, total_heal.data, LineWidth=2)
ylabel('$\char"0023$ of cases', 'Interpreter', 'latex')
% title('\textbf{Observed vs Total Healed Population}', 'Interpreter', 'latex')
grid on
% lgd = legend('Estimated Total', 'Observed','Interpreter', 'latex','location','southeast');
% lgd.FontSize = 18;
yax = ylabel('$\mathbf{\char"0023}$ of cases','Interpreter','latex');
% yax.FontSize = 14;
set(gca, 'TickLabelInterpreter', 'Latex')
xlim([dataAUS.Date(1), dataAUS.Date(end)])
box on


% "I" Undetected infected, estimated model

var_area.x = [dataAUS.Date; flip(dataAUS.Date)];
var_area.y = [dark_plus_avg.data'; flip(dark_minus_avg.data')];
var_area.Ivar = [infected_avg_plus; flip(infected_avg_minus)];

figure(13)
plot(dataAUS.Date, dark_avg.data, LineWidth=2, Color=[0.8500 0.3250 0.0980])
hold on
fill(var_area.x,var_area.y,[0.8500 0.3250 0.0980],'FaceAlpha',.3,'EdgeColor','none')
hold on
yax = ylabel('$\mathbf{\char"0023}$ of cases','Interpreter','latex');
% yax.FontSize = 14;
% title('\textbf{Estimated Hidden Cases per day}','Interpreter','latex')
grid on
lgd = legend('Hidden Cases','95 \% Confidence Interval','Interpreter','latex','Location','northeast');
% lgd.FontSize = 18;
xlim([dataAUS.Date(1), dataAUS.Date(end)])
set(gca, 'TickLabelInterpreter', 'Latex')
box on

save('Var_InfectedAUS.mat', 'var_area');

% Ratio between "ESTIMATED TOTAL INFECTED / DETECTED INFECTED" - Trend

figure(14)
plot(dataAUS.Date, ratio_undetected, LineWidth=2, Color=[0.8500 0.3250 0.0980])
% title('\textbf{Ratio $\zeta$ Undetected - Detected Infections }','Interpreter','latex')
grid on
% legend('$\zeta$','Interpreter','latex')
yax = ylabel('Ratio Level','Interpreter','latex');
% yax.FontSize = 14;
xlim([dataAUS.Date(1), dataAUS.Date(end)])
set(gca, 'TickLabelInterpreter', 'Latex')


% Actual I group

figure(20)
plot(dataAUS.Date, DataStructAUS.Filterd.ActiveCases./ratio_undetected2', LineWidth=1.5, LineStyle="--",Color=[0 0 0])
hold on
plot(dataAUS.Date, infected_avg, LineWidth=1.5)
hold on
fill(var_area.x,var_area.Ivar,[0.8500 0.3250 0.0980],'FaceAlpha',.2,'EdgeColor','none')
ylabel('$\char"0023$ of cases','Interpreter','latex')
grid on
legend('Real Data - Diagnosed', 'Dark Data - Undetected','95\% Confidence Interval','Interpreter','latex')
xlim([dataAUS.Date(1), dataAUS.Date(end)])
set(gca, 'TickLabelInterpreter', 'Latex')

%% Saving all data

DataStructAUS.Filterd.Vax = Vax_Data.people_fully_vaccinated;

SIDTTHE_AUS = {infected_avg, 'Infected'; DataStructAUS.Filterd.ActiveCases, 'Diagnosed'; DataStructAUS.Filterd.Hospitalised, 'Hospitalised'; DataStructAUS.Filterd.ICUs, 'ICU threatend'; DataStructAUS.Filterd.Deaths, 'Deceased'; total_heal.data, 'HealedAug'; DataStructAUS.Filterd.Recoveries, 'Healed'};
filename = 'SIDTTHE_AUS.mat';
save(filename, 'SIDTTHE_AUS');

filename2 = 'DataStructAUS.mat';
save(filename2, 'DataStructAUS');

%% Australian Variants

variants = readtable('/Users/marcodelloro/Desktop/Thesis/MSc-Thesis-TUe/Data_Collection/Australian Dataset/variantsAUS 2.xlsx');

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
data_barplot = removevars(variants,'Var1');
data_barplot(:, 11:18) = []; % in oder to have less variants (not interesting)
data_barplot = table2array(data_barplot)*100;
date_barplot = variants.Var1;

figure()
h = bar(date_barplot, data_barplot, 0.90, "stacked");
for ii = 1:size(customColors, 1)
    set(h(1, ii), 'facecolor', customColors{ii, 1})
end
ylabel('$\%$ Variant', 'Interpreter', 'latex')
legend('B.1.1.7 - ALPHA VARIANT', 'B.1.351 - BETA VARIANT', 'P.1 - GAMMA VARIANT', 'B.1.617.2 - DELTA VARIANT', ...
       'BA.1 - OMICRON BA.1', 'BA.2 - OMICRON BA.2', 'BA.2.12.1', 'BA.4 - OMICRON BA.4', 'BA.5 - OMICRON BA.5', 'BQ.1', ...
       'BA.2.86', 'Others',...
       'Interpreter', 'latex', 'Location', 'southoutside', 'NumColumns', 3)
xlim([variants.Var1(1)-5, variants.Var1(end)+5])
ylim([0, 100])
set(gca, 'TickLabelInterpreter', 'Latex')

% Smoothing of Most importat variants to have a "Variants Plot"

xii = [0:1:29];
p=1;

SmoothVar.date = variants.Var1;
SmoothVar.alpha = csaps(xii,variants.("delta"),p,xii);
SmoothVar.omi_BA1 = csaps(xii,variants.("omicron_BA_1_"),p,xii);
SmoothVar.omi_BA2 = csaps(xii,variants.("omicron_BA_2_"),p,xii);
SmoothVar.omi_BA4 = csaps(xii,variants.("omicron_BA_4_"),p,xii);
SmoothVar.omi_BA5 = csaps(xii,variants.("Omicron_BA_5_"),p,xii);

save('SmoothVariantasAus.mat','SmoothVar');
