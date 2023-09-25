% The following code prepares data gathered in an external code to be fitted in the
% SIDTTHE (simplified SIDARTHE) MODEL.
% The goal is to have data describing the current situation, day by day (no cumulative data)

clear all
close all
clc

load("cleandata.mat")

% Modified SIDARTHE model - SIDTTHE 
% S = Susceptible Healthy People
% I = Infected people, not yet detected (Untested)
% D = Detected, both symptomatic and asymptomatic
% T1 = Currently hospitalised individuals
% T2 = Currently IC individuals
% H = Currently Healed individuals
% E = Currently Expired individuals 

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


