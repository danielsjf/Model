clc;
clear all;

year = '2013';

%% Spot prices
file_f = ['spotmarket_data_',year,'.xls'];
[num,txt,raw] = xlsread(file_f);
spot = num(:,1:24);
spot = reshape(spot, 1, numel(spot));

save(['Spot',year,'.mat'], 'spot');