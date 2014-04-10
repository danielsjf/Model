%% postprocessing illustration
clc
clear all
close all

load('WPFE_data.mat'); 

TSwind.forecast =forecast;
TSwind.production = measurement;
TSwind.total_cap = total_cap;
TSwind.delta_T = 0.25;
TSwind.n_samples = length(forecast);


day = 1; 
load(['day',num2str(day),'.mat']);
scen_error = scen_error';
[samples_scenario,n_scen_gen] = size(scen_error);
forecast = TSwind.forecast((day-1)*samples_scenario+1:day*samples_scenario);
measurement = TSwind.production((day-1)*samples_scenario+1:day*samples_scenario);

scen_error_f = zeros(samples_scenario,n_scen_gen);

for a = 1:100
   scen_error_f(:,a) = scen_error(:,a)+forecast; 
end

time = [0.25:0.25:24]';

data = [time forecast measurement scen_error_f];
save('illustration_scen_gen.dat','data','-ascii');

load('WPFE_ee.mat')
load('WPFE_ep.mat')
[scen_prob_mat,scen_prob] = scen_probability(scen_error,n_scen_gen,forecast,pdf_stable,samples_scenario,pdf_stable_ee,w_int,w_ee,'calc');

