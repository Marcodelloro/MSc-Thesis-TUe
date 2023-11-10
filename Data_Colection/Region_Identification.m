%% Different Pandemics Region Representation

% The following code is for the identification of the zones based on the different data from the pandemics;
% The plot is representing a 2-axis plot of the Diagnosed and the ICUs recovered, which are (as asssumed) 
% the two most important parameters to understand the effect of the pandemic.

clc
clear all
close all

load('/Users/marcodelloro/Desktop/Thesis/MSc-Thesis-TUe/Data_Colection/SIDTTHE_data.mat');

Npop = 59240329; % Total Population of Italy
ICUs = 4000; % Total number of ICUs in Italy @ 9 Ottobre 2020

D_data = SIDTTHE_data{2,1}.data / Npop;     % D group - Diagnosed Population
T2_data = SIDTTHE_data{4,1}.data /ICUs ;    % T2 group - ICUs (Severly Threatned) Population
        
N = 399; % if we are using the DAILY time discretization, then N = 399
      
% New time vector (daily data) - daily discretization
new_time = linspace(SIDTTHE_data{1, 1}.date(1), SIDTTHE_data{1, 1}.date(end), N); % 399 data points (daily data)

DataArray = {D_data, T2_data};

% Upsampling of the dataset using Interp1, so that we can handle better the divided day by day
for ii = 1:numel(DataArray)

    original_data = DataArray{1, ii};
    upsampled_DataArray{ii,1} = interp1(SIDTTHE_data{1, 1}.date, original_data, new_time);

end

%% Plot of the behaviour of the states evolving in the different areas of the 

perc_ICUs = [0.1 0.3 0.65 1];

area.yellow.x = [0 0 perc_ICUs(1) perc_ICUs(1)];
area.yellow.y = [0 max(upsampled_DataArray{1,1})*1.5 max(upsampled_DataArray{1,1})*1.5 0];

area.orange.x = [perc_ICUs(1) perc_ICUs(1) perc_ICUs(2) perc_ICUs(2)];
area.orange.y = [0 max(upsampled_DataArray{1,1})*1.25 max(upsampled_DataArray{1,1})*1.25 0];

area.red.x = [perc_ICUs(2) perc_ICUs(2) perc_ICUs(3) perc_ICUs(3)];
area.red.y = [0 max(upsampled_DataArray{1,1})*1.05 max(upsampled_DataArray{1,1})*1.05 0];

area.red_max.x = [perc_ICUs(3) perc_ICUs(3) perc_ICUs(4) perc_ICUs(4)];
area.red_max.y = [0 max(upsampled_DataArray{1,1})*1.15 max(upsampled_DataArray{1,1})*1.15 0];

% Additional Areas

area.red_max2.x = [perc_ICUs(2) perc_ICUs(2) perc_ICUs(3) perc_ICUs(3)];
area.red_max2.y = [max(upsampled_DataArray{1,1})*1.05 max(upsampled_DataArray{1,1})*1.15 max(upsampled_DataArray{1,1})*1.15 max(upsampled_DataArray{1,1})*1.05];

area.red_max3.x = [perc_ICUs(1) perc_ICUs(1) perc_ICUs(2) perc_ICUs(2)];
area.red_max3.y = [max(upsampled_DataArray{1,1})*1.25 max(upsampled_DataArray{1,1})*1.5 max(upsampled_DataArray{1,1})*1.5 max(upsampled_DataArray{1,1})*1.25];

area.red_max4.x = [0 0 1 1 0.3 0.3];
area.red_max4.y = [max(upsampled_DataArray{1,1})*1.5 0.02 0.02 max(upsampled_DataArray{1,1})*1.15 max(upsampled_DataArray{1,1})*1.15 max(upsampled_DataArray{1,1})*1.5 ];


customColorsAreas = {   [1, 0.6, 0.2],
                        [1, 0.3, 0.05],
                        [1, 0.2, 0.05],
                        [0.8, 0.2, 0.1] 
                    };

% Plot of the behaviour of the states evolving

policy_idx = [1 40 69 116 141 190 242 294 399];

customColors2 = {   [1, 0.6, 0.2],
                    [1, 0.3, 0.05],
                    [1, 0.2, 0.05],
                    [0.8, 0.2, 0.1],
                    [1, 0.2, 0.05],            
                    [0.8, 0.2, 0.1],
                    [1, 0.3, 0.05],
                    [1, 0.6, 0.2]
                };

markers = { 'o', 's', 'h', 'd', 'v', 'd', 'h', 's', 'o'};

figure()
for jj=1:length(perc_ICUs)-1
    xline(perc_ICUs(jj),":",'HandleVisibility', 'off')
    hold on 
end
fill(area.yellow.x, area.yellow.y, customColorsAreas{1,1} ,'FaceAlpha',.5,'EdgeColor', 'none','HandleVisibility', 'off')
hold on
fill(area.orange.x, area.orange.y, customColorsAreas{2,1} ,'FaceAlpha',.5,'EdgeColor', 'none','HandleVisibility', 'off')
hold on
fill(area.red.x, area.red.y, customColorsAreas{3,1} ,'FaceAlpha',.5,'EdgeColor', 'none','HandleVisibility', 'off')
hold on 
fill(area.red_max2.x, area.red_max2.y, customColorsAreas{4,1} ,'FaceAlpha',.5,'EdgeColor', 'none','HandleVisibility', 'off')
hold on 
fill(area.red_max3.x, area.red_max3.y, customColorsAreas{4,1} ,'FaceAlpha',.5,'EdgeColor', 'none','HandleVisibility', 'off')
hold on 
fill(area.red_max4.x, area.red_max4.y, customColorsAreas{4,1} ,'FaceAlpha',.5,'EdgeColor', 'none','HandleVisibility', 'off')
hold on
fill(area.red_max.x, area.red_max.y, customColorsAreas{4,1} ,'FaceAlpha',.5,'EdgeColor', 'none','HandleVisibility', 'off')
for ii = 1:length(policy_idx)-1
    scatter( upsampled_DataArray{2,1}(policy_idx(ii):policy_idx(ii+1)), upsampled_DataArray{1,1}(policy_idx(ii):policy_idx(ii+1)), "filled", 'MarkerFaceColor', customColors2{ii},'MarkerEdgeColor','k')
    hold on
end
ylabel('Percentage of Diagnosed','Interpreter','latex')
xlabel('Intensive Care Units Occupied','Interpreter','latex')
title('\textbf{\textit{D} and pressure on ICus - Comparison}','Interpreter','latex')
grid on
ylim([0 0.02])
set(gca, 'TickLabelInterpreter', 'Latex')





