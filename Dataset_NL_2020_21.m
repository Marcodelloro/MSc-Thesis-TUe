clc
clear all 
close all

% The data gathered https://data.rivm.nl/meta/srv/dut/catalog.search#/search?topicCat=health
% are spanning from june 2020 to october 2021; to our interest are the 57
% weeks between 31 august 2020 and 3 october 2021, so the year of
% COVID-19 outbrake, taking into account 2021 second wave and intial
% vaccination of the population of Netherlands as January 2021.

% here are gathered all sorts of data from RIVM, classified and polotted.
% Only few of them are actually used to fit the whole SIDTTHE (SIMPLIFIED
% SIDARTHE) model, the rest are to have a complete overview of the problem
% at the moment

%% Total tests + Positive tested 

data_test = readtable('/Users/marcodelloro/Desktop/Thesis/Codes/Netherlands_dataset.xlsx','Sheet','COVID-19_Testing_Positive');
data_ICU= readtable('/Users/marcodelloro/Desktop/Thesis/Codes/Netherlands_dataset.xlsx','Sheet','COVID-19_ICUs_age');
data_ICU = data_ICU(487:end,:);
data_ReD = readtable('/Users/marcodelloro/Desktop/Thesis/Codes/Netherlands_dataset.xlsx','Sheet','COVID-19_RepCases_Deaths'); % reported cases and Deceased
data_ReD = data_ReD(68821:end,:);
dataset_vax= readtable('/Users/marcodelloro/Desktop/Thesis/Codes/Netherlands_dataset.xlsx','Sheet','COVID-19_vaccination');
dataset_variants = readtable('/Users/marcodelloro/Desktop/Thesis/Codes/Netherlands_dataset.xlsx','Sheet','Variants','Range','A:S');

% loading of the measures coming from the government
policy = readtable('/Users/marcodelloro/Desktop/Thesis/Codes/Netherlands_dataset.xlsx','Sheet','GovActions','Range','B25:B37');
policy.DateOfActions = datetime(policy.DateOfActions, 'InputFormat', 'dd MMMM yyyy');

%% Total tests + Positive tested 
data_test.Date_of_statistics = datetime(data_test.Date_of_statistics, 'InputFormat', 'dd-MMM-yyyy HH:mm:ss');
total_tests = groupsummary(data_test, "Date_of_statistics", 'sum', 'Tested_with_result');
positive_tests = groupsummary(data_test, "Date_of_statistics", 'sum', 'Tested_positive');

total_tests = total_tests(92:end,:);  % Only gather tests from monday 31 August
positive_tests = positive_tests(92:end,:);

total_tests= removevars(total_tests, 'GroupCount');
positive_tests= removevars(positive_tests, 'GroupCount');
window = 7; % weekly average

for ii = 1:window:size(total_tests,1)
    tt_avg.data(ii)= mean(total_tests.sum_Tested_with_result(ii:ii + window -1));
    pt_avg.data(ii) = mean(positive_tests.sum_Tested_positive(ii:ii+ window-1));
end

tt_avg.data = nonzeros(tt_avg.data);
pt_avg.data = nonzeros(pt_avg.data);

tt_avg.csum = cumsum(tt_avg.data,1);
pt_avg.csum = cumsum(pt_avg.data,1);

% Creation of new data time in weeks

t1 = datetime(2020,8,31);
t2 = datetime(2021,10,3);
tt_avg.date = t1:caldays(7):t2;
pt_avg.date = t1:caldays(7):t2;

% As Kohler does, a normalization over whole population must be brought on
Npop = 17475415; %Total population of Netherlands in 2021

% Not normalized data figures
figure()
plot(tt_avg.date, tt_avg.data, 'LineWidth', 1.5)
hold on 
plot(pt_avg.date, pt_avg.data, 'LineWidth', 1.5)
ylabel('$\char"0023$ of cases', 'Interpreter', 'latex')
title('\textbf{Total tests vs Positive tests per week}', 'Interpreter', 'latex')
grid on
legend('Total', 'Positive','Interpreter', 'latex')
xlim([tt_avg.date(1), tt_avg.date(end)])
set(gca, 'TickLabelInterpreter', 'Latex')

figure()
plot(tt_avg.date, tt_avg.csum, 'LineWidth',1.5)
hold on 
plot(pt_avg.date, pt_avg.csum, 'LineWidth',1.5)
ylabel('$\char"0023$ of cases', 'Interpreter', 'latex')
title('\textbf{Cumulate Total tests vs Cumulate Positive tests per week}','Interpreter', 'latex')
grid on
legend('Total', 'Positive','Interpreter', 'latex','Location','northwest')
xlim([tt_avg.date(1), tt_avg.date(end)])
set(gca, 'TickLabelInterpreter', 'Latex')

% As Kohler does, a normalization over whole population must be brought on
Npop = 17475415; %Total population of Netherlands in 2021

tt_avg.data_norm = tt_avg.data/Npop;
pt_avg.data_norm = pt_avg.data/Npop;

tt_avg.csum_norm = cumsum(tt_avg.data_norm,1);
pt_avg.csum_norm = cumsum(pt_avg.data_norm,1);

% Normalized data figures
figure()
plot(tt_avg.date, tt_avg.data_norm, 'LineWidth',1.5)
hold on 
plot(pt_avg.date, pt_avg.data_norm, 'LineWidth',1.5)
ylabel('$\%$ of cases over Population', 'Interpreter', 'latex')
title('\textbf{$\%$ Total tests vs $\%$ Positive tests per week}','Interpreter', 'latex')
grid on
legend('Total', 'Positive','Interpreter', 'latex')
xlim([tt_avg.date(1), tt_avg.date(end)])
set(gca, 'TickLabelInterpreter', 'Latex')

% Save only useful data
cleandata.testing.total.date = tt_avg.date;
cleandata.testing.total.data = tt_avg.data_norm;

cleandata.testing.positive.date = pt_avg.date;
cleandata.testing.positive.data = pt_avg.data_norm;

figure()
plot(tt_avg.date, tt_avg.csum_norm, 'LineWidth',1.5)
hold on 
plot(pt_avg.date, pt_avg.csum_norm, 'LineWidth',1.5)
ylabel('$\%$ of cases over Population', 'Interpreter', 'latex')
title('\textbf{Cumulate $\%$ Tests vs Cumulate $\%$ Positive tests per week}','Interpreter', 'latex')
grid on
legend('Total', 'Positive','Interpreter', 'latex','Location','northwest')
xlim([tt_avg.date(1), tt_avg.date(end)])
set(gca, 'TickLabelInterpreter', 'Latex')



%% ICU ADMISSIONS

% Division based on age
child = data_ICU(strcmp(data_ICU.Age_group, '0-14'), :);
teen = data_ICU(strcmp(data_ICU.Age_group, '15-19'), :);
adults = data_ICU(ismember(data_ICU.Age_group, {'20-24','30-34', '35-39', '40-44'}), :);
midlife = data_ICU(ismember(data_ICU.Age_group, {'50-54', '50-54', '55-59'}), :);
senior = data_ICU(ismember(data_ICU.Age_group, {'60-64', '65-69', '70-74'}), :);
elderly = data_ICU(ismember(data_ICU.Age_group, {'75-79', '80-84', '85-89'}), :);
sengen = data_ICU(ismember(data_ICU.Age_group, {'90+'}), :);
unkwn = data_ICU(ismember(data_ICU.Age_group, {'Unknown'}), :);

% Sum all present weeks
child_summed.t1 = groupsummary(child, 'Date_of_statistics_week_start', 'sum', 'Hospital_admission');
child_summed.t2 = groupsummary(child, 'Date_of_statistics_week_start', 'sum', 'IC_admission');
teen_summed.t1 = groupsummary(teen, 'Date_of_statistics_week_start', 'sum', 'Hospital_admission');
teen_summed.t2 = groupsummary(teen, 'Date_of_statistics_week_start', 'sum', 'IC_admission');
adults_summed.t1 = groupsummary(adults, 'Date_of_statistics_week_start', 'sum', 'Hospital_admission');
adults_summed.t2 = groupsummary(adults, 'Date_of_statistics_week_start', 'sum', 'IC_admission');
midlife_summed.t1 = groupsummary(midlife, 'Date_of_statistics_week_start', 'sum', 'Hospital_admission');
midlife_summed.t2 = groupsummary(midlife, 'Date_of_statistics_week_start', 'sum', 'IC_admission');
senior_summed.t1 = groupsummary(senior, 'Date_of_statistics_week_start', 'sum', 'Hospital_admission');
senior_summed.t2 = groupsummary(senior, 'Date_of_statistics_week_start', 'sum', 'IC_admission');
elderly_summed.t1 = groupsummary(elderly, 'Date_of_statistics_week_start', 'sum', 'Hospital_admission');
elderly_summed.t2 = groupsummary(elderly, 'Date_of_statistics_week_start', 'sum', 'IC_admission');
sengen_summed.t1 = groupsummary(sengen, 'Date_of_statistics_week_start', 'sum', 'Hospital_admission');
sengen_summed.t2 = groupsummary(sengen, 'Date_of_statistics_week_start', 'sum', 'IC_admission');
unkwn_summed.t1 = groupsummary(unkwn, 'Date_of_statistics_week_start', 'sum', 'Hospital_admission');
unkwn_summed.t2 = groupsummary(unkwn, 'Date_of_statistics_week_start', 'sum', 'IC_admission');

% Cumulated sum of all the age groups t1 and t2
child_summed.t1.csum = cumsum(child_summed.t1.sum_Hospital_admission,1);
child_summed.t2.csum = cumsum(child_summed.t2.sum_IC_admission,1);
teen_summed.t1.csum = cumsum(teen_summed.t1.sum_Hospital_admission,1);
teen_summed.t2.csum = cumsum(teen_summed.t2.sum_IC_admission,1);
adults_summed.t1.csum = cumsum(adults_summed.t1.sum_Hospital_admission,1);
adults_summed.t2.csum = cumsum(adults_summed.t2.sum_IC_admission,1);
midlife_summed.t1.csum = cumsum(midlife_summed.t1.sum_Hospital_admission,1);
midlife_summed.t2.csum = cumsum(midlife_summed.t2.sum_IC_admission,1);
senior_summed.t1.csum = cumsum(senior_summed.t1.sum_Hospital_admission,1);
senior_summed.t2.csum = cumsum(senior_summed.t2.sum_IC_admission,1);
elderly_summed.t1.csum = cumsum(elderly_summed.t1.sum_Hospital_admission,1);
elderly_summed.t2.csum = cumsum(elderly_summed.t2.sum_IC_admission,1);
sengen_summed.t1.csum = cumsum(elderly_summed.t1.sum_Hospital_admission,1);
sengen_summed.t2.csum = cumsum(elderly_summed.t2.sum_IC_admission,1);
unkwn_summed.t1.csum = cumsum(elderly_summed.t1.sum_Hospital_admission,1);
unkwn_summed.t2.csum = cumsum(elderly_summed.t2.sum_IC_admission,1);

% For the sake of the fitting operations all data differentiated by age
% group will be summed to have a normalization over the whole population

ICUs_Total = child_summed.t2.sum_IC_admission + teen_summed.t2.sum_IC_admission + adults_summed.t2.sum_IC_admission + ...
             midlife_summed.t2.sum_IC_admission + senior_summed.t2.sum_IC_admission + elderly_summed.t2.sum_IC_admission + ...
             sengen_summed.t2.sum_IC_admission + unkwn_summed.t2.sum_IC_admission;
Hosp_Total = child_summed.t1.sum_Hospital_admission + teen_summed.t1.sum_Hospital_admission + adults_summed.t1.sum_Hospital_admission + ...
             midlife_summed.t1.sum_Hospital_admission + senior_summed.t1.sum_Hospital_admission + elderly_summed.t1.sum_Hospital_admission + ...
             sengen_summed.t1.sum_Hospital_admission + unkwn_summed.t1.sum_Hospital_admission;
ICUs_Total = ICUs_Total/Npop; % Normalization over population
Hosp_Total = Hosp_Total/Npop;

%saving of the data in 'cleandata' struct
cleandata.IC_Hos.ICUs.date = child_summed.t1.Date_of_statistics_week_start;
cleandata.IC_Hos.ICUs.data = ICUs_Total;
cleandata.IC_Hos.Hosp.date = child_summed.t1.Date_of_statistics_week_start;
cleandata.IC_Hos.Hosp.data = Hosp_Total;

% Not normalized data figures - Trend plot 
figure()
plot(child_summed.t1.Date_of_statistics_week_start, child_summed.t1.sum_Hospital_admission, 'LineWidth',1.5)
hold on
plot(child_summed.t1.Date_of_statistics_week_start, teen_summed.t1.sum_Hospital_admission, 'LineWidth',1.5)
hold on
plot(child_summed.t1.Date_of_statistics_week_start, adults_summed.t1.sum_Hospital_admission, 'LineWidth',1.5)
hold on
plot(child_summed.t1.Date_of_statistics_week_start, midlife_summed.t1.sum_Hospital_admission, 'LineWidth',1.5)
hold on
plot(child_summed.t1.Date_of_statistics_week_start, senior_summed.t1.sum_Hospital_admission, 'LineWidth',1.5)
hold on
plot(child_summed.t1.Date_of_statistics_week_start, elderly_summed.t1.sum_Hospital_admission, 'LineWidth',1.5)
hold on
plot(child_summed.t1.Date_of_statistics_week_start, sengen_summed.t1.sum_Hospital_admission, 'LineWidth',1.5)
hold on
plot(child_summed.t1.Date_of_statistics_week_start, unkwn_summed.t1.sum_Hospital_admission, 'LineWidth',1.5)
ylabel('$\char"0023$ of cases', 'Interpreter', 'latex')
title('\textbf{Hospitalisation Trends Based on Age Groups}', 'Interpreter', 'latex')
grid on
legend('Child', 'Teen', 'Adult', 'Midlife', 'Senior', 'Elderly', 'Senoctogenarian', 'Unknown','Interpreter', 'latex')
xlim([tt_avg.date(1), tt_avg.date(end)])
set(gca, 'TickLabelInterpreter', 'Latex')


figure()
plot(child_summed.t1.Date_of_statistics_week_start, child_summed.t2.sum_IC_admission, 'LineWidth',1.5)
hold on
plot(child_summed.t1.Date_of_statistics_week_start, teen_summed.t2.sum_IC_admission, 'LineWidth',1.5)
hold on
plot(child_summed.t1.Date_of_statistics_week_start, adults_summed.t2.sum_IC_admission, 'LineWidth',1.5)
hold on
plot(child_summed.t1.Date_of_statistics_week_start, midlife_summed.t2.sum_IC_admission, 'LineWidth',1.5)
hold on
plot(child_summed.t1.Date_of_statistics_week_start, senior_summed.t2.sum_IC_admission, 'LineWidth',1.5)
hold on
plot(child_summed.t1.Date_of_statistics_week_start, elderly_summed.t2.sum_IC_admission, 'LineWidth',1.5)
hold on
plot(child_summed.t1.Date_of_statistics_week_start, sengen_summed.t2.sum_IC_admission, 'LineWidth',1.5)
hold on
plot(child_summed.t1.Date_of_statistics_week_start, unkwn_summed.t2.sum_IC_admission, 'LineWidth',1.5)
ylabel('$\char"0023$ of cases', 'Interpreter', 'latex')
title('\textbf{ICU Trends Based on Age Groups}', 'Interpreter', 'latex')
grid on
legend('Child', 'Teen', 'Adult', 'Midlife', 'Senior', 'Elderly', 'Senoctogenarian', 'Unknown','Interpreter', 'latex')
set(gca, 'TickLabelInterpreter', 'Latex')
xlim([tt_avg.date(1), tt_avg.date(end)])



% Not normalized data figures - Cumulated plot 

figure()
plot(child_summed.t1.Date_of_statistics_week_start, child_summed.t1.csum, 'LineWidth',1.5)
hold on
plot(child_summed.t1.Date_of_statistics_week_start, teen_summed.t1.csum, 'LineWidth',1.5)
hold on
plot(child_summed.t1.Date_of_statistics_week_start, adults_summed.t1.csum, 'LineWidth',1.5)
hold on
plot(child_summed.t1.Date_of_statistics_week_start, midlife_summed.t1.csum, 'LineWidth',1.5)
hold on
plot(child_summed.t1.Date_of_statistics_week_start, senior_summed.t1.csum, 'LineWidth',1.5)
hold on
plot(child_summed.t1.Date_of_statistics_week_start, elderly_summed.t1.csum, 'LineWidth',1.5)
hold on
plot(child_summed.t1.Date_of_statistics_week_start, sengen_summed.t1.csum, 'LineWidth',1.5)
hold on
plot(child_summed.t1.Date_of_statistics_week_start, unkwn_summed.t1.csum, 'LineWidth',1.5)
ylabel('$\char"0023$ of cases', 'Interpreter', 'latex')
title('\textbf{Hospitalisation Cumulated Data Based on Age Groups}', 'Interpreter', 'latex')
grid on
legend('Child', 'Teen', 'Adult', 'Midlife', 'Senior', 'Elderly', 'Senoctogenarian', 'Unknown','Interpreter', 'latex','Location','northeastoutside')
xlim([tt_avg.date(1), tt_avg.date(end)])
set(gca, 'TickLabelInterpreter', 'Latex')


figure()
plot(child_summed.t1.Date_of_statistics_week_start, child_summed.t2.csum, 'LineWidth',1.5)
hold on
plot(child_summed.t1.Date_of_statistics_week_start, teen_summed.t2.csum, 'LineWidth',1.5)
hold on
plot(child_summed.t1.Date_of_statistics_week_start, adults_summed.t2.csum, 'LineWidth',1.5)
hold on
plot(child_summed.t1.Date_of_statistics_week_start, midlife_summed.t2.csum, 'LineWidth',1.5)
hold on
plot(child_summed.t1.Date_of_statistics_week_start, senior_summed.t2.csum, 'LineWidth',1.5)
hold on
plot(child_summed.t1.Date_of_statistics_week_start, elderly_summed.t2.csum, 'LineWidth',1.5)
hold on
plot(child_summed.t1.Date_of_statistics_week_start, sengen_summed.t2.csum, 'LineWidth',1.5)
hold on
plot(child_summed.t1.Date_of_statistics_week_start, unkwn_summed.t2.csum, 'LineWidth',1.5)
ylabel('$\char"0023$ of cases', 'Interpreter', 'latex')
title('\textbf{ICUs Cumulated Data Based on Age Groups}', 'Interpreter', 'latex')
grid on
legend('Child', 'Teen', 'Adult', 'Midlife', 'Senior', 'Elderly', 'Senoctogenarian', 'Unknown','Interpreter', 'latex','Location','northeastoutside')
xlim([tt_avg.date(1), tt_avg.date(end)])
set(gca, 'TickLabelInterpreter', 'Latex')


%% Total reported + Reported Deaths

clear  rep_avg.data dec_avg.data heal_avg.data rep_avg.csum dec_avg.csum heal_avg.csum

data_ReD.Date_of_publication = datetime(data_ReD.Date_of_publication, 'InputFormat', 'dd-MMM-yyyy HH:mm:ss');
total_rep = groupsummary(data_ReD, "Date_of_publication", 'sum', 'Total_reported'); % Total reported cases
total_dec = groupsummary(data_ReD, "Date_of_publication", 'sum', 'Deceased');  % Total reported deceased

total_rep= removevars(total_rep, 'GroupCount');
total_dec= removevars(total_dec, 'GroupCount');
total_heal.Total_healed = total_rep.sum_Total_reported - total_dec.sum_Deceased;
total_heal.Date_of_publication = total_rep.Date_of_publication;
window = 7; % weekly average

for ii = 1:window:size(total_rep,1)
    rep_avg.data(ii) = mean(total_rep.sum_Total_reported(ii:ii + window -1));
    dec_avg.data(ii) = mean(total_dec.sum_Deceased(ii:ii+ window-1));
    heal_avg.data(ii) = mean(total_heal.Total_healed(ii:ii+ window-1));
end

rep_avg.data = nonzeros(rep_avg.data);
dec_avg.data = nonzeros(dec_avg.data);
heal_avg.data = nonzeros(heal_avg.data);

rep_avg.csum = cumsum(rep_avg.data,1);
dec_avg.csum = cumsum(dec_avg.data,1);
heal_avg.csum = cumsum(heal_avg.data,1);

% Creation of new data time in weeks

t1 = datetime(2020,8,31);
t2 = datetime(2021,10,3);
rep_avg.date = t1:caldays(7):t2;
dec_avg.date = t1:caldays(7):t2;
heal_avg.date = t1:caldays(7):t2;

% Comparison polt to see if data gathered before are coherent TOTAL POSITIVE vs TOTAL REPORTED
figure()
plot(pt_avg.date, pt_avg.data, 'LineWidth',1.5)
hold on
plot(rep_avg.date, rep_avg.data, 'LineWidth',1.5)
ylabel('$\char"0023$ of cases','Interpreter','latex')
title('\textbf{Comparison between Different Data Sources}', 'Interpreter','latex')
grid on
legend('Positive from tests', 'Reported from Data', 'Interpreter','latex')
xlim([tt_avg.date(1), tt_avg.date(end)])
set(gca, 'TickLabelInterpreter', 'Latex')

figure()
plot(dec_avg.date, dec_avg.data, 'LineWidth',1.5)
title('\textbf{Total Deceased per week}', 'Interpreter','latex')
grid on
legend('Deceased Cases', 'Interpreter','latex')
xlim([tt_avg.date(1), tt_avg.date(end)])
set(gca, 'TickLabelInterpreter', 'Latex')

figure()
plot(dec_avg.date, dec_avg.csum, 'LineWidth',1.5)
ylabel('$\char"0023$ of cases','Interpreter','latex')
title('\textbf{Cumulate Total Deceased per week}','Interpreter','latex')
grid on
legend('Deceased Cases','Interpreter','latex','Location','southeast')
xlim([tt_avg.date(1), tt_avg.date(end)])
set(gca, 'TickLabelInterpreter', 'Latex')

figure()
plot(heal_avg.date, heal_avg.data, 'LineWidth',1.5)
title('\textbf{Total Healed per week}','Interpreter','latex')
grid on
legend('Healed Cases','Interpreter','latex')
xlim([tt_avg.date(1), tt_avg.date(end)])
set(gca, 'TickLabelInterpreter', 'Latex')

figure()
plot(heal_avg.date, heal_avg.csum, 'LineWidth',1.5)
ylabel('$\char"0023$ of cases','Interpreter','latex')
title('\textbf{Cumulate Total Healed per week}','Interpreter','latex')
grid on
legend('Healed Cases','Interpreter','latex','Location','southeast')
xlim([tt_avg.date(1), tt_avg.date(end)])
set(gca, 'TickLabelInterpreter', 'Latex')

% Normalise data over the whole population

rep_avg.data_norm = rep_avg.data/Npop;
dec_avg.data_norm = dec_avg.data/Npop;
heal_avg.data_norm = heal_avg.data/Npop;

rep_avg.csum_norm = cumsum(rep_avg.data_norm,1);
dec_avg.csum_norm = cumsum(dec_avg.data_norm,1);
heal_avg.csum_norm = cumsum(heal_avg.data_norm,1);

% Normalized data figures
figure()
plot(rep_avg.date, rep_avg.data_norm, 'LineWidth',1.5)
ylabel('$\%$ of cases over Population','Interpreter','latex')
title('\textbf{$\%$ Reported Positive Cases per week}','Interpreter','latex')
grid on
legend('Reported Positive','Interpreter','latex')
xlim([tt_avg.date(1), tt_avg.date(end)])
set(gca, 'TickLabelInterpreter', 'Latex')

figure()
plot(dec_avg.date, dec_avg.data_norm, 'LineWidth',1.5)
ylabel('$\%$ of cases over Population','Interpreter','latex')
title('\textbf{$\%$ Deceased Cases per week}','Interpreter','latex')
grid on
legend('Deceased','Interpreter','latex')
xlim([tt_avg.date(1), tt_avg.date(end)])
set(gca, 'TickLabelInterpreter', 'Latex')

figure()
plot(heal_avg.date, heal_avg.data_norm, 'LineWidth',1.5)
ylabel('$\%$ of cases over Population','Interpreter','latex')
title('\textbf{$\%$ Healed Cases per week}','Interpreter','latex')
grid on
legend('Healed','Interpreter','latex')
xlim([tt_avg.date(1), tt_avg.date(end)])
set(gca, 'TickLabelInterpreter', 'Latex')

% Saving of the normalized data in 'cleandata' struct
cleandata.Reported.Pos.date = rep_avg.date;
cleandata.Reported.Pos.data = rep_avg.data_norm;
cleandata.Reported.Dec.date = dec_avg.date;
cleandata.Reported.Dec.data = dec_avg.data_norm;
cleandata.Reported.Heal.date = heal_avg.date;
cleandata.Reported.Heal.data = heal_avg.data_norm;

figure()
plot(rep_avg.date, rep_avg.csum_norm, 'LineWidth',1.5)
title('\textbf{Cumulate $\%$ Reported Positive Cases per week}','Interpreter','latex')
grid on
legend('Reported Positive','Interpreter','latex','Location','southeast')
xlim([tt_avg.date(1), tt_avg.date(end)])
set(gca, 'TickLabelInterpreter', 'Latex')
ylabel('$\%$ of cases over Population','Interpreter','latex')

figure()
plot(dec_avg.date, dec_avg.csum_norm, 'LineWidth',1.5)
ylabel('$\%$ of cases over Population','Interpreter','latex')
title('\textbf{Cumulate $\%$ Deceased Cases per week}','Interpreter','latex')
grid on
legend('Deceased','Interpreter','latex','Location','southeast')
xlim([tt_avg.date(1), tt_avg.date(end)])
set(gca, 'TickLabelInterpreter', 'Latex')

figure()
plot(heal_avg.date, heal_avg.csum_norm, 'LineWidth',1.5)
ylabel('$\%$ of cases over Population','Interpreter','latex')
title('\textbf{Cumulate $\%$ Healed Cases per week}','Interpreter','latex')
grid on
legend('Healed','Interpreter','latex','Location','southeast')
xlim([tt_avg.date(1), tt_avg.date(end)])
set(gca, 'TickLabelInterpreter', 'Latex')


%% Vaccination activities 

dataset_vax.Coverage_primary_partly = str2double(dataset_vax.Coverage_primary_partly);
dataset_vax.Coverage_primary_completed = str2double(dataset_vax.Coverage_primary_completed);
dataset_vax.Coverage_repeat_vaccination_autumn_round = str2double(dataset_vax.Coverage_repeat_vaccination_autumn_round);

% Division based on birth age
under12 = dataset_vax(strcmp(dataset_vax.Birth_year, '<=2011'), :);
over60 = dataset_vax(strcmp(dataset_vax.Birth_year, '<=1963'), :);
over18 = dataset_vax(strcmp(dataset_vax.Birth_year, '<2006'), :);
teenagers = dataset_vax(strcmp(dataset_vax.Birth_year, '2006-2011'), :);

for ii=1:size(under12,1)
    
    under12.tot_dose1_partial(ii,:) = round(under12.Populatie(ii,:) * under12.Coverage_primary_partly(ii,:)/100);
    over60.tot_dose1_partial(ii,:) = round(over60.Populatie(ii,:) * over60.Coverage_primary_partly(ii,:)/100);
    teenagers.tot_dose1_partial(ii,:) = round(teenagers.Populatie(ii,:) * teenagers.Coverage_primary_partly(ii,:)/100);
    over18.tot_dose1_partial(ii,:) = round(over18.Populatie(ii,:) * over18.Coverage_primary_partly(ii,:)/100);

    under12.tot_dose1_full(ii,:) = round(under12.Populatie(ii,:) * under12.Coverage_primary_completed(ii,:)/100);
    over60.tot_dose1_full(ii,:) = round(over60.Populatie(ii,:) * over60.Coverage_primary_completed(ii,:)/100);
    teenagers.tot_dose1_full(ii,:) = round(teenagers.Populatie(ii,:) * teenagers.Coverage_primary_completed(ii,:)/100);
    over18.tot_dose1_full(ii,:) = round(over18.Populatie(ii,:) * over18.Coverage_primary_completed(ii,:)/100);

    under12.tot_booster(ii,:) = round(under12.Populatie(ii,:) * under12.Coverage_repeat_vaccination_autumn_round(ii,:)/100);
    over60.tot_booster(ii,:) = round(over60.Populatie(ii,:) * over60.Coverage_repeat_vaccination_autumn_round(ii,:)/100);
    teenagers.tot_booster(ii,:) = round(teenagers.Populatie(ii,:) * teenagers.Coverage_repeat_vaccination_autumn_round(ii,:)/100);
    over18.tot_booster(ii,:) = round(over18.Populatie(ii,:) * over18.Coverage_repeat_vaccination_autumn_round(ii,:)/100);

end

% Now sum (still by agegroup) all the days to have a total of the whole Netherlands

under12_summed.partial = groupsummary(under12, "Date_of_statistics", 'sum', 'tot_dose1_partial');
under12_summed.full = groupsummary(under12, "Date_of_statistics", 'sum', 'tot_dose1_full');
under12_summed.booster = groupsummary(under12, "Date_of_statistics", 'sum', 'tot_booster');
under12_summed.pop = groupsummary(under12, "Date_of_statistics", 'sum', 'Populatie');

over60_summed.partial = groupsummary(over60, "Date_of_statistics", 'sum', 'tot_dose1_partial');
over60_summed.full = groupsummary(over60, "Date_of_statistics", 'sum', 'tot_dose1_full');
over60_summed.booster = groupsummary(over60, "Date_of_statistics", 'sum', 'tot_booster');
over60_summed.pop = groupsummary(over60, "Date_of_statistics", 'sum', 'Populatie');

teenagers_summed.partial = groupsummary(teenagers, "Date_of_statistics", 'sum', 'tot_dose1_partial');
teenagers_summed.full = groupsummary(teenagers, "Date_of_statistics", 'sum', 'tot_dose1_full');
teenagers_summed.booster = groupsummary(teenagers, "Date_of_statistics", 'sum', 'tot_booster');
teenagers_summed.pop = groupsummary(teenagers, "Date_of_statistics", 'sum', 'Populatie');

over18_summed.partial = groupsummary(over18, "Date_of_statistics", 'sum', 'tot_dose1_partial');
over18_summed.full = groupsummary(over18, "Date_of_statistics", 'sum', 'tot_dose1_full');
over18_summed.booster = groupsummary(over18, "Date_of_statistics", 'sum', 'tot_booster');
over18_summed.pop = groupsummary(over18, "Date_of_statistics", 'sum', 'Populatie');


% Plots as a mattter of vaccines per group age - PARTIAL COVER

figure()
plot(under12_summed.partial.Date_of_statistics, under12_summed.partial.sum_tot_dose1_partial, LineWidth=1.5)
hold on
plot(teenagers_summed.partial.Date_of_statistics, teenagers_summed.partial.sum_tot_dose1_partial, LineWidth=1.5)
hold on
plot(over18_summed.partial.Date_of_statistics, over18_summed.partial.sum_tot_dose1_partial, LineWidth=1.5)
hold on
plot(over60_summed.partial.Date_of_statistics, over60_summed.partial.sum_tot_dose1_partial, LineWidth=1.5)
grid on
ylabel('$\char"0023$ of doses','Interpreter','latex')
legend('Child', 'Teen', 'Adult', 'Senior','Interpreter','latex','Location','northeastoutside')
title('\textbf{Total Number of Vaccines Allocated per Age Group - Partial Cover}','Interpreter','latex')
set(gca, 'TickLabelInterpreter', 'Latex')
xlim([over60_summed.partial.Date_of_statistics(1), over60_summed.partial.Date_of_statistics(end)])

figure()
plot(under12_summed.partial.Date_of_statistics, under12_summed.partial.sum_tot_dose1_partial./under12_summed.pop.sum_Populatie, LineWidth=1.5)
hold on
plot(teenagers_summed.partial.Date_of_statistics, teenagers_summed.partial.sum_tot_dose1_partial./teenagers_summed.pop.sum_Populatie, LineWidth=1.5)
hold on
plot(over18_summed.partial.Date_of_statistics, over18_summed.partial.sum_tot_dose1_partial./over18_summed.pop.sum_Populatie, LineWidth=1.5)
hold on
plot(over60_summed.partial.Date_of_statistics, over60_summed.partial.sum_tot_dose1_partial./over60_summed.pop.sum_Populatie, LineWidth=1.5)
grid on
ylabel('$\%$ of Population','Interpreter','latex')
legend('Child', 'Teen', 'Adult', 'Senior','Interpreter','latex','Location','northeastoutside')
title('\textbf{Percentage of Vaccines Allocated per Age Group - Partial Cover}','Interpreter','latex')
set(gca, 'TickLabelInterpreter', 'Latex')
xlim([over60_summed.partial.Date_of_statistics(1), over60_summed.partial.Date_of_statistics(end)])


% Plots as a mattter of vaccines per group age - TOTAL COVER

figure()
plot(under12_summed.full.Date_of_statistics, under12_summed.full.sum_tot_dose1_full, LineWidth=1.5)
hold on
plot(teenagers_summed.full.Date_of_statistics, teenagers_summed.full.sum_tot_dose1_full, LineWidth=1.5)
hold on
plot(over18_summed.full.Date_of_statistics, over18_summed.full.sum_tot_dose1_full, LineWidth=1.5)
hold on
plot(over60_summed.full.Date_of_statistics, over60_summed.full.sum_tot_dose1_full, LineWidth=1.5)
grid on
ylabel('$\char"0023$ of doses','Interpreter','latex')
legend('Child', 'Teen', 'Adult', 'Senior','Interpreter','latex','Location','northeastoutside')
title('\textbf{Total Number of Vaccines Allocated per Age Group - Total Cover}','Interpreter','latex')
set(gca, 'TickLabelInterpreter', 'Latex')
xlim([over60_summed.partial.Date_of_statistics(1), over60_summed.partial.Date_of_statistics(end)])


figure()
plot(under12_summed.full.Date_of_statistics, under12_summed.full.sum_tot_dose1_full./under12_summed.pop.sum_Populatie, LineWidth=1.5)
hold on
plot(teenagers_summed.full.Date_of_statistics, teenagers_summed.full.sum_tot_dose1_full./teenagers_summed.pop.sum_Populatie, LineWidth=1.5)
hold on
plot(over18_summed.full.Date_of_statistics, over18_summed.full.sum_tot_dose1_full./over18_summed.pop.sum_Populatie, LineWidth=1.5)
hold on
plot(over60_summed.full.Date_of_statistics, over60_summed.full.sum_tot_dose1_full./over60_summed.pop.sum_Populatie, LineWidth=1.5)
grid on
ylabel('$\%$ of Population','Interpreter','latex')
legend('Child', 'Teen', 'Adult', 'Senior','Interpreter','latex','Location','northeastoutside')
title('\textbf{Percentage of Vaccines Allocated per Age Group - Total Cover}','Interpreter','latex')
set(gca, 'TickLabelInterpreter', 'Latex')
xlim([over60_summed.partial.Date_of_statistics(1), over60_summed.partial.Date_of_statistics(end)])


% Plots as a mattter of vaccines per group age - BOOSTER COVER

figure()
plot(under12_summed.booster.Date_of_statistics, under12_summed.booster.sum_tot_booster, LineWidth=1.5)
hold on
plot(teenagers_summed.booster.Date_of_statistics, teenagers_summed.booster.sum_tot_booster, LineWidth=1.5)
hold on
plot(over18_summed.booster.Date_of_statistics, over18_summed.booster.sum_tot_booster, LineWidth=1.5)
hold on
plot(over60_summed.booster.Date_of_statistics, over60_summed.booster.sum_tot_booster, LineWidth=1.5)
grid on
ylabel('$\char"0023$ of doses','Interpreter','latex')
legend('Child', 'Teen', 'Adult', 'Senior','Interpreter','latex','Location','northeastoutside')
title('\textbf{Total Number of Vaccines Allocated per Age Group - Booster Dose}','Interpreter','latex')
set(gca, 'TickLabelInterpreter', 'Latex')
xlim([over60_summed.partial.Date_of_statistics(1), over60_summed.partial.Date_of_statistics(end)])


figure()
plot(under12_summed.booster.Date_of_statistics, under12_summed.booster.sum_tot_booster./under12_summed.pop.sum_Populatie, LineWidth=1.5)
hold on
plot(teenagers_summed.booster.Date_of_statistics, teenagers_summed.booster.sum_tot_booster./teenagers_summed.pop.sum_Populatie, LineWidth=1.5)
hold on
plot(over18_summed.full.Date_of_statistics, over18_summed.full.sum_tot_dose1_full./over18_summed.pop.sum_Populatie, LineWidth=1.5)
hold on
plot(over18_summed.booster.Date_of_statistics, over18_summed.booster.sum_tot_booster./over60_summed.pop.sum_Populatie, LineWidth=1.5)
grid on
ylabel('$\%$ of Population','Interpreter','latex')
legend('Child', 'Teen', 'Adult', 'Senior','Interpreter','latex','Location','northeastoutside')
title('\textbf{Percentage of Vaccines Allocated per Age Group - Booster Dose}','Interpreter','latex')
set(gca, 'TickLabelInterpreter', 'Latex')
xlim([over60_summed.partial.Date_of_statistics(1), over60_summed.partial.Date_of_statistics(end)])


%% SARS-CoV-V2 Variants

for ii = 1:size(dataset_variants,1)

    % Estraction year from the string
    data_str = char(dataset_variants.Category);
    year_str = data_str(ii,1:4);
    year = str2double(year_str);
    
    % Estraction week from the string
    week_str = data_str(ii,7:end);
    week = str2double(week_str);

    newdate = datetime(year, 1, 1) + calweeks(week - 1) + days(6);
    dataset_variants.newdate(ii) = newdate;
end

% Figures 
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
                [0 1 0],    
                [1 0 1],
                [0 1 1],
                [0.7, 0.3, 0.5],  
                [0.3, 0.6, 0.9],   
                [0.4, 0.8, 0.2],     
                [0.7, 0.7, 0.7]};

figure()
for ii = 2:19
    area(dataset_variants.newdate, dataset_variants{:,ii},'FaceColor', [customColors{ii-1,1}], 'FaceAlpha', 0.7);
    hold on
end

grid on 
title('\textbf{SARS-CoV-2 Variants in The Netherlands}','Interpreter','latex')
legend('SARS-CoV-2', 'B.1.1.7 - ALPHA VARIANT', 'B.1.351', 'P.1', 'B.1.617.2 - DELTA VARIANT', 'BA.1 - OMICRON BA1 VARIANT', 'BA.2 - OMICRON BA2 VARIANT', 'BA.2.12.1','BA.4','BA.5 - OMICRON BA5 VARIANT','BQ.1 - OMICRON BQ1 VARIANT	','BA.2.75','XBB','XBB.1.5 OMICRON XBB VARIANT','XBB.1.9',  'XBF','XBB.1.16','EG.5','Interpreter','latex','Location','northeastoutside')
xlim([dataset_variants.newdate(1), dataset_variants.newdate(end)])
set(gca, 'TickLabelInterpreter', 'Latex')
ylabel('$\%$ Variant','Interpreter','latex')


% modifiy data for the plot
data_barplot = removevars(dataset_variants,'Category');
data_barplot = removevars(data_barplot,'newdate');
data_barplot = table2array(data_barplot);
date_barplot = dataset_variants{:,20};

figure()
h = bar(date_barplot,data_barplot, 3, "stacked");
for ii = 1:size(customColors, 1)
    h(1, ii).FaceColor = [customColors{ii,1}];
end

ylabel('$\%$ Variant','Interpreter','latex')
title('\textbf{SARS-CoV-2 Variants in The Netherlands - Trend per Week}','Interpreter','latex')
legend('SARS-CoV-2', 'B.1.1.7 - ALPHA VARIANT', 'B.1.351', 'P.1', 'B.1.617.2 - DELTA VARIANT', 'BA.1 - OMICRON BA1 VARIANT', 'BA.2 - OMICRON BA2 VARIANT', 'BA.2.12.1','BA.4','BA.5 - OMICRON BA5 VARIANT','BQ.1 - OMICRON BQ1 VARIANT','BA.2.75','XBB','XBB.1.5 OMICRON XBB VARIANT','XBB.1.9','XBF','XBB.1.16','EG.5','Interpreter','latex','Location','northeastoutside')
xlim([dataset_variants.newdate(1), dataset_variants.newdate(end)])
ylim([0,100.1])
set(gca, 'TickLabelInterpreter', 'Latex')

%%  Adding pictures with areas highlighting the Action periods
% 11 different areas will be covered

policy_dates = vertcat(tt_avg.date(1), policy.DateOfActions);

customColors2 = {   [1, 0.9, 0.6],
                    [1, 0.8, 0.4],
                    [1, 0.6, 0.2],
                    [1, 0.4, 0.1],
                    [1, 0.3, 0.05],
                    [1, 0.2, 0.05],
                    [0.8, 0.2, 0.1],
                    [1, 0.2, 0.05],            
                    [1, 0.3, 0.05],              
                    [1, 0.4, 0.1],
                    [1, 0.6, 0.2],
                    [1, 0.8, 0.4]
                    
                };

for ii = 1:size(policy_dates)-1
    area.x(ii, :) = [policy_dates(ii) policy_dates(ii) policy_dates(ii+1) policy_dates(ii+1)];
    area.y(ii, :) = [0 7e-4 7e-4 0];
    area.y2(ii,:) = [0 7e-6 7e-6 0];
end

figure()
for ii = 1:size(policy_dates)-1
    fill(area.x(ii, :) ,area.y(1, :), customColors2{ii,1} ,'FaceAlpha',.5,'EdgeColor', 'none','HandleVisibility', 'off')
    hold on
    xline(policy_dates(ii),":",'HandleVisibility', 'off')
    hold on 
end
hold on
plot(rep_avg.date, rep_avg.data_norm, 'k','LineWidth',1.5, 'DisplayName', 'Reported Positive')
ylabel('$\%$ of cases over Population','Interpreter','latex')
title('\textbf{$\%$ Reported Positive Cases per week}','Interpreter','latex')
grid on
legend('Interpreter','latex')
xlim([tt_avg.date(1), tt_avg.date(end)])
ylim([0, 7e-4])
set(gca, 'TickLabelInterpreter', 'Latex')


figure()
for ii = 1:size(policy_dates)-1
    fill(area.x(ii, :) ,area.y2(1, :), customColors2{ii,1} ,'FaceAlpha',.5,'EdgeColor', 'none','HandleVisibility', 'off')
    hold on
    xline(policy_dates(ii),":",'HandleVisibility', 'off')
    hold on 
end
hold on
plot(dec_avg.date, dec_avg.data_norm,'k','LineWidth',1.5, 'DisplayName', 'Deceased')
ylabel('$\%$ of cases over Population','Interpreter','latex')
title('\textbf{$\%$ Deceased Cases per week}','Interpreter','latex')
grid on
legend('Interpreter','latex')
xlim([tt_avg.date(1), tt_avg.date(end)])
ylim([0, 7e-6])
set(gca, 'TickLabelInterpreter', 'Latex')

close all

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

% H(t) =  [ delta_N(t) * ( delta_N(t) -1 ) ] / [ 1 + delta_N(t-1) - delta_D(t) ] 
% H(t) = Hidden cases bias-corrected form by Chao


% NOTE: EVEN IF THE PAPER PRODUCES HIS DATA ON DAILY BASIS, I STICK WITH WEEKS, SO DATA CAN BE MORE SMOOTH;
% MAYBE IN THIS CASE WOULD BE BETETR TO ADD THE HEALED GROUP ASWELL IN THE 'COUNTED TWICE CASES', BUT CAN BE FIXED IN THE FUTURE

delta_N = rep_avg.data; 
delta_D = dec_avg.data;

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
    disp (ii)
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

% Plots of the variance and the mean value

var_area.x = [rep_avg.date flip(rep_avg.date)];
var_area.y = [H_plus flip(H_minus)];

figure()
fill(var_area.x,var_area.y,[0 0.4470 0.7410],'FaceAlpha',.5,'EdgeColor', 'none','HandleVisibility', 'off')
hold on
plot(rep_avg.date, H, LineWidth=1.5, Color=[0 0.4470 0.7410])
ylabel('$\char"0023$ of cases','Interpreter','latex')
title('\textbf{Estimated Hidden Cases per day}','Interpreter','latex')
grid on
legend('Hidden Cases','Interpreter','latex','Location','northeast')
xlim([tt_avg.date(1), tt_avg.date(end)])
set(gca, 'TickLabelInterpreter', 'Latex')

ratio_observed = (H' + rep_avg.data) ./ rep_avg.data ;

figure()
plot(rep_avg.date, ratio_observed, LineWidth=1.5)
ylabel('Ratio Total to Observed','Interpreter','latex')
title('\textbf{Ratio Total to Observed Trend}','Interpreter','latex')
grid on
legend('Ratio T/O','Interpreter','latex','Location','southeast')
xlim([tt_avg.date(1), tt_avg.date(end)])
set(gca, 'TickLabelInterpreter', 'Latex')


%% Estimate of the total number of Healed Population (H)

% According to CSSEGISandData (John Hopkins University Covid-19 World Database) Netherlands
% does not publish recovered data.


%% Gathering of the data foe every compartment
% 
% filename = 'cleandata.mat';
% save(filename, '-struct', 'cleandata'); 