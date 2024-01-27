%% DATA FROM ITALY
% Vaccine dataset from the Italian dataset https://github.com/italia/covid19-opendata-vaccini/tree/master
% Our interest are the 57 weeks between 31 august 2020 and 3 october 2021, so the year of
% COVID-19 outbrake, taking into account 2021 second wave and intial vaccination of the population.

clc
clear all 
close all


% Note that for the time window used, only 1st dose and 2nd dose are relevant
vaxRaw20 = readtable('somministrazioni-vaccini-latest-2020.csv');
vaxRaw21 = readtable('somministrazioni-vaccini-latest-2021.csv');

vaxRaw20.data = datetime(vaxRaw20.data, 'InputFormat', 'yyyy-MM-dd');
vaxRaw21.data = datetime(vaxRaw21.data, 'InputFormat', 'yyyy-MM-dd');

selectedColumns = {'data','d1', 'd2', 'dpi', 'db1'};

vaxRaw20 = vaxRaw20(:, selectedColumns);
vaxRaw21 = vaxRaw21(:, selectedColumns);

vax20 = groupsummary(vaxRaw20, 'data', 'sum');
vax21 = groupsummary(vaxRaw21, 'data', 'sum');

vax_tot = [table2array(vax20(:,3:end)); table2array(vax21(1:276,3:end))];
missing_days = 399 - size(vax_tot,1);

wait_period = zeros(missing_days,4);
% wait_period = ones(missing_days,4) * 0.1;

vax_tot = [wait_period; vax_tot];

% recreating the final table so can be manged better

t1 = datetime(2020,8,31);
t2 = datetime(2021,10,3);

date = t1:caldays(1):t2;

columnNames = {'d1', 'd2', 'dpi', 'db1'};
vaxData = array2table(vax_tot, 'VariableNames', columnNames);

vaxData.date = date';
vaxData = vaxData(:, ['date', columnNames]);

% Actual cumulated sum of double-dosed
vaxData.sum_d2 = cumsum(vaxData.d2);
vaxData.sum_dpi = cumsum(vaxData.dpi);

% Vaccines Cumulated 
figure()
plot(vaxData.date, vaxData.sum_d2, LineWidth=1.5)
hold on
plot(vaxData.date, vaxData.sum_dpi, LineWidth=1.5)
ylabel('$\char"0023$ of cases','Interpreter','latex')
title('\textbf{Vaccine Compartment Trend}','Interpreter','latex')
grid on
legend('Vaccines 2nd dose','Post Infection Dose','latex', 'Location','northwest')
xlim([date(1), date(end)])
set(gca, 'TickLabelInterpreter', 'Latex')

% Vaccines Day-by-Day
figure()
plot(vaxData.date, vaxData.d1, LineWidth=1.5)
hold on
plot(vaxData.date, vaxData.d2, LineWidth=1.5)
hold on 
plot(vaxData.date, vaxData.dpi, LineWidth=1.5)
hold on 
plot(vaxData.date, vaxData.db1, LineWidth=1.5)
ylabel('$\char"0023$ of cases','Interpreter','latex')
title('\textbf{Raw Vaccine Data}','Interpreter','latex')
grid on
legend('1st Dose', '2nd Dose', 'Infected Single Dose', '1st Booster', '2nd Booster','latex', 'Location','northeast')
xlim([date(1), date(end)])
set(gca, 'TickLabelInterpreter', 'Latex')


%% Smoothing and filtering
p = 0.4;
xi = [0:1:398];
ys = csaps(xi,vaxData.d1,p,xi);

% Apply a Savitzky - Golay filter on the cubic data
order = 3;
framelen = 17;
sgf = sgolayfilt(ys,order,framelen);

figure()
plot(vaxData.date, vaxData.d1, LineWidth=1.5, LineStyle="--")
hold on
plot(xi, ys, LineWidth=2)
hold on
plot(xi, sgf, LineWidth=2)
grid on
legend('Raw', 'Cubic Spline', 'Cubic Spline + SGF','latex', 'Location','northeast')
xlim([date(1), date(end)])
ylim([0, 5.5e5])
set(gca, 'TickLabelInterpreter', 'Latex')

%% Saving in .mat

filename = 'vaxData.mat';
save(filename, 'vaxData'); 

filename2 = 'Filtered VaxData';
