% postprocessing 

clc
clear all
close all

load('WPFE_data.mat'); 
load('WPFE_ep.mat');
load('WPFE_ee.mat');

TSwind.forecast =forecast;
TSwind.production = measurement;
TSwind.total_cap = total_cap;
TSwind.delta_T = 0.25;
TSwind.n_samples = length(forecast);
% % Prediction errors | Forecast
n_int = 1/w_int;
ll_int = [0:w_int:1-w_int];
ul_int = [0+w_int:w_int:1];
ll_error = [-1:w_int:1-w_int];
ul_error = [-1+w_int:w_int:1];
Bs_tot = zeros(365,11,7);
vec_val = [-0.5:0.1:0.5];

for day = 1:365
disp(['day ',num2str(day)])
load(['day',num2str(day),'.mat']);
%scen_error = scen_error';
[n_scen_gen,samples_scenario] = size(scen_error);
hist_error = measurement((day-1)*96+1:day*96)-forecast((day-1)*96+1:day*96);
Bs = zeros(samples_scenario,11,7);

% functional: values
count = 1;
for d = 0:4:24
if d == 0
    d = 1;
end
for c = 1:11
for b = 1:samples_scenario
temp = zeros(n_scen_gen,1);
for a = 1:n_scen_gen
    temp(a) = functional_var(scen_error(a,:),b,d,vec_val(c));    
end
Pt = 1/n_scen_gen*sum(temp);
obs = functional_var(hist_error,b,d,vec_val(c));    
Bs(b,c,count) = (Pt-obs)^2;
end
Bs_tot(day,c,count) = 1/(samples_scenario)*sum(Bs(:,c,count));
end
count = count+1;
end
end

Bs_final = zeros(11,7);
Bs_final_var = zeros(11,7);
for c = 1:11
    for count = 1:7
    Bs_final(c,count) = 1/365*sum(Bs_tot(:,c,count));
    Bs_final_var(c,count) = var(Bs_tot(:,c,count));
    end
end

