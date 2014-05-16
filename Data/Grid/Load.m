clc;
clear all;

year = '2013';

%% Spot prices
file_f = ['Load_Forecast_',year,'.xls'];
[num,txt,raw] = xlsread(file_f,1);
LoadFor = num(1:end-2,1:24);

file_f = ['Generation_Forecast_Historical_By_Units_',year,'.xls'];
[num,txt,raw] = xlsread(file_f,1);
LoadForUnit = num(3:end,49);
LoadForUnit = [LoadForUnit; LoadForUnit(end); LoadForUnit(end)]; % no data for 30 and 31 december 2013
LoadForUnit(isnan(LoadForUnit)) = 0;
LoadForUnit = repmat(LoadForUnit,1,24);

file_f = ['Installed_Power_Historical_',year,'.xls'];
[num,txt,raw] = xlsread(file_f,1);
InstallCap = num(:,1);
InstallCap = [InstallCap; InstallCap(end)]; % no data for 31 december 2013
InstallCap = repmat(InstallCap,1,24);

DifCap = InstallCap - LoadFor;

LoadFor = reshape(LoadFor', 1, numel(LoadFor));
LoadForUnit = reshape(LoadForUnit', 1, numel(LoadForUnit));
InstallCap = reshape(InstallCap', 1, numel(InstallCap));
DifCap = reshape(DifCap', 1, numel(DifCap));

figure()
plot(LoadFor)
figure()
plot(LoadForUnit)
figure()
plot(InstallCap)
figure()
plot(DifCap)
% spot(isnan(spot))=mean(spot(isfinite(spot)));
% 
save(['Load',year,'.mat'], 'LoadFor', 'InstallCap', 'DifCap', 'LoadForUnit');