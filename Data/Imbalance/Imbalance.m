clc;
clear all;

year = '2013';

%% Imbalance
imbalance = [];
i = 1;
while i <= 12
    file_f = ['Imbalance-',year,'-',sprintf('%02d', i),'.xls'];
    [num,txt,raw] = xlsread(file_f);
    imbalance = [imbalance; num];
    i = i+1;
end

NRV = imbalance(:,1);
SI = imbalance(:,2);
alpha = imbalance(:,3);
MIP = imbalance(:,4);
MDP = imbalance(:,5);
POS = imbalance(:,6);
NEG = imbalance(:,7);

save(['Imbalance',year,'.mat'], 'NRV', 'SI', 'alpha', 'MIP', 'MDP', 'POS', 'NEG');