%% DATA FROM ITALY
% Vaccine dataset from the Italian dataset https://github.com/italia/covid19-opendata-vaccini/tree/master
% Our interest are the 57 weeks between 31 august 2020 and 3 october 2021, so the year of
% COVID-19 outbrake, taking into account 2021 second wave and intial vaccination of the population.

clc
clear all 
close all

vaxRaw20 = readtable('somministrazioni-vaccini-latest-2020.csv');
vaxRaw21 = readtable('somministrazioni-vaccini-latest-2021.csv');

vaxRaw20.data = datetime(vaxRaw20.data, 'InputFormat', 'yyyy-MM-dd');
vaxRaw21.data = datetime(vaxRaw21.data, 'InputFormat', 'yyyy-MM-dd');

selectedColumns = {'data','d1', 'd2', 'dpi', 'db1', 'db2', 'db3'};

vaxRaw20 = vaxRaw20(:, selectedColumns);
vaxRaw21 = vaxRaw21(:, selectedColumns);

vax20 = groupsummary(vaxRaw20, 'data', 'sum');
vax21 = groupsummary(vaxRaw21, 'data', 'sum');

vax_tot = [table2array(vax20(:,3:end)); table2array(vax21(1:276,3:end))];
missing_days = 399 - size(vax_tot,1);

wait_period = zeros(missing_days,6);

vax_tot = [wait_period; vax_tot];

% recreating the final table so can be manged better

t1 = datetime(2020,8,31);
t2 = datetime(2021,10,3);

date = t1:caldays(1):t2;

columnNames = {'d1', 'd2', 'dpi', 'db1', 'db2', 'db3'};
vaxData = array2table(vax_tot, 'VariableNames', columnNames);

vaxData.date = date';
vaxData = vaxData(:, ['date', columnNames]);

%% Saving in .mat

filename = 'vaxData.mat';
save(filename, 'vaxData'); 
