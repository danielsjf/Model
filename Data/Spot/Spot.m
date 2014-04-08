clc;
clear all;

year = '2013';

%% Spot prices
file_f = ['spotmarket_data_',year,'.xls'];
[num,txt,raw] = xlsread(file_f,2);
spot = num(:,1:24);
spot = reshape(spot, 1, numel(spot));
spot(isnan(spot))=mean(spot(isfinite(spot)));

save(['Spot',year,'.mat'], 'spot');