clc
clear all 
close all

% The data gathered https://data.rivm.nl/meta/srv/dut/catalog.search#/search?topicCat=health
% are spanning from june 2020 to october 2021; to our interest are the 57
% weeks between 31 august 2020 and 3 october 2021, so the year of
% COVID-19 outbrake, taking into account 2021 second wave and intial
% vaccination of the population of Netherlands as January 2021.

%% Total tests + Positive tested 

data_test = readtable('2020_21_Covid19_NL.xlsx','Sheet','COVID-19_Testing_Positive');

data_test.Date_of_statistics = datetime(data_test.Date_of_statistics, 'InputFormat', 'dd-MMM-yyyy HH:mm:ss');
total_tests = groupsummary(data_test, "Date_of_statistics", 'sum', 'Tested_with_result');
positive_tests = groupsummary(data_test, "Date_of_statistics", 'sum', 'Tested_positive');

total_tests = total_tests(92:end,:);  % Only gather tests from monday 31 August
positive_tests = positive_tests(92:end,:);

total_tests= removevars(total_tests, 'GroupCount');
positive_tests= removevars(positive_tests, 'GroupCount');
window = 7; % weekly average

for ii = 1:window:size(total_tests,1)
    disp(ii)
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
plot(tt_avg.date, tt_avg.data, 'LineWidth',1.5)
hold on 
plot(pt_avg.date, pt_avg.data, 'LineWidth',1.5)
ylabel('# of cases')
title('Total tests vs Positive tests per week')
grid on
legend('Total', 'Positive')
xlim([tt_avg.date(1), tt_avg.date(end)])

figure()
plot(tt_avg.date, tt_avg.csum, 'LineWidth',1.5)
hold on 
plot(pt_avg.date, pt_avg.csum, 'LineWidth',1.5)
ylabel('# of cases')
title('Cumulate Total tests vs Cumulate Positive tests per week')
grid on
legend('Total', 'Positive')
xlim([tt_avg.date(1), tt_avg.date(end)])

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
ylabel('% of cases over Population')
title('% Total tests vs % Positive tests per week')
grid on
legend('Total', 'Positive')
xlim([tt_avg.date(1), tt_avg.date(end)])

figure()
plot(tt_avg.date, tt_avg.csum_norm, 'LineWidth',1.5)
hold on 
plot(pt_avg.date, pt_avg.csum_norm, 'LineWidth',1.5)
ylabel('% of cases over Population')
title('Cumulate % Tests vs Cumulate % Positive tests per week')
grid on
legend('Total', 'Positive')
xlim([tt_avg.date(1), tt_avg.date(end)])


%% ICU ADMISSIONS

data_ICU= readtable('2020_21_Covid19_NL.xlsx','Sheet','COVID-19_ICUs_age');
data_ICU = data_ICU(487:end,:);

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


% Not normalized data figures
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
ylabel('# of cases')
title('Hospitalisation Trends Based on Age Groups')
grid on
legend('Child', 'Teen', 'Adults', 'Midlife', 'Senior', 'Elderly', 'Senoctogenarian', 'Unknown' )
xlim([tt_avg.date(1), tt_avg.date(end)])


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
ylabel('# of cases')
title('ICU Trends Based on Age Groups')
grid on
legend('Child', 'Teen', 'Adults', 'Midlife', 'Senior', 'Elderly', 'Senoctogenarian', 'Unknown' )
xlim([tt_avg.date(1), tt_avg.date(end)])

