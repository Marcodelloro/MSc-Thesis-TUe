function [DataArray, new_time, matrixData] = LoadingData(SIDTTHE_data,vaxData)

Npop = 59240329; % Total Population of Italy

new_time = SIDTTHE_data{1,1}.date;

I_data = SIDTTHE_data{1,1}.data / Npop;                             % I group - Infected Population **
D_data = SIDTTHE_data{2,1}.data / Npop;                             % D group - Diagnosed Population
T1_data = SIDTTHE_data{3,1}.data / Npop;                            % T1 group - Hospitalised (Lightly Threatned) Population
T2_data = SIDTTHE_data{4,1}.data / Npop;                            % T2 group - ICUs (Severly Threatned) Population
H_data = (SIDTTHE_data{6,1}.data - vaxData.sum_dpi') / Npop;        % T2 group - Healed Population **
E_data = SIDTTHE_data{5,1}.data  / Npop;                            % E group - Expired (Deceased) Population
V_data = (vaxData.sum_d2 + vaxData.sum_dpi)' / Npop;
S_data = ones(length(I_data),1)' - (I_data + D_data + T1_data + T2_data + H_data + E_data + V_data ); 

% ** the groups I and H are based on strong assumptions and dark numbers estimates

DataArray = {S_data ; I_data ; D_data ; T1_data ; T2_data ; H_data ; E_data ; V_data };
matrixData = [S_data ; I_data ; D_data ; T1_data ; T2_data ; H_data ; E_data ; V_data];

end