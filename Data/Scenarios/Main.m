%% Scenario generation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This script uses the input from the excel sheet to generate different
% scenarios for different days. The output is a matlab file with the
% original scenarios and their probabilities and the reduced scenarios and
% their probabilities.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Some parameters for the simulations
% start & end day for the simulation
start_day = xlsread('simulation_parameters.xlsx','Sheet1','C12'); % First day of the simulation
end_day = xlsread('simulation_parameters.xlsx','Sheet1','C13'); % Last day of the simulation

% TP
time_period = xlsread('simulation_parameters.xlsx','Sheet1','C31'); % Quarters per hour (=0.25)

% Load load profile
load(['INPUT' filesep 'LOAD.mat'],'demand_net','other'); % Read out electricity demand (net and other)

% Load the historical time series TSWind.mat
load(['INPUT' filesep 'WPFE_data.mat']);  % Read out wind turbine data (TSwind: production and forecast) (WPFE = wind power forecast error)

% the annual share of wind
share_wind = xlsread('simulation_parameters.xlsx','Sheet1','C14');
% rescale wind + scenarios
tot_wind = sum(TSwind.production);
tot_load = sum(sum(demand_net+other));
factor_wind = share_wind/(tot_wind/tot_load);
cap_wind = factor_wind*TSwind.tot_cap;

% select the relevant days
TSwind.forecast = TSwind.forecast((start_day-1)*96+1:end_day*96);
TSwind.production = TSwind.production((start_day-1)*96+1:end_day*96);
demand_net = demand_net(start_day:end_day,:);
other = other(start_day:end_day,:);

%number of days
n_days = end_day-start_day+1;

for day = 1:n_days
if day > 1
disp('          ')
end
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp(['%%%%%%%%%%%%%%%%%%%     Day ',num2str(day),'    %%%%%%%%%%%%%%%%%%%%%%%'])
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('          ')

%% Generate some scenarios - for now just for testing: load results
% forecast scenario
if day == 1
scen_forecast = TSwind.forecast(1+ 96*(day-1):96*day);
else
scen_forecast = [TSwind.forecast(1+96*(day-2):96*(day-1))+actual_realization; TSwind.forecast(1+ 96*(day-1):96*day)];    
end
% actual realization
actual_realization = TSwind.production(1+96*(day-1):96*day)-TSwind.forecast(1+ 96*(day-1):96*day);

% Load the statistical analysis of the WPFE. The stable distribution is used to generate the probabilities of the various scenarios.
% Probability of an error given the forecast
load(['INPUT' filesep 'WPFE_ep.mat']);
cdf_stable = cumsum(pdf_stable')';
% Probability of an error given the previous error
load(['INPUT' filesep 'WPFE_ee.mat']);
% length_scenario: length of the scenarios 
length_scenario = xlsread('simulation_parameters.xlsx','Sheet1','C20');
samples_scenario = length_scenario*1/TSwind.delta_T;
% Generation 
n_scen_gen = xlsread('simulation_parameters.xlsx','Sheet1','C21');
MU = zeros(1,samples_scenario);
coef_fit_var = xlsread('simulation_parameters.xlsx','Sheet1','C23:Q23');
norm_var = xlsread('simulation_parameters.xlsx','Sheet1','C24');
% Reduction
n_scen_select = 30;
[useless,ff_mode] = xlsread('simulation_parameters.xlsx','Sheet1','C25');

% Scenario generaton and reduction
if day == 1
[scen_lim,scen_red] = SCEN(n_scen_gen,MU,coef_fit_var,norm_var,scen_forecast,samples_scenario,cdf_stable,pdf_stable,pdf_stable_ee,w_int,w_ee,error_vec,n_scen_select-1,actual_realization,ff_mode,day);
else
[scen_lim,scen_red] = SCEN(n_scen_gen,MU,coef_fit_var,norm_var,scen_forecast(97:192),samples_scenario,cdf_stable,pdf_stable,pdf_stable_ee,w_int,w_ee,error_vec,n_scen_select-1,actual_realization,ff_mode,day);    
end

save(['SCENARIOS' filesep 'day',num2str(day)],'scen_lim','scen_red');

end
