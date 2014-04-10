%% INPUT selecting days   
% needed for ieee 

clear all
close all
clc

time_period = 0.25; 

load(['INPUT' filesep 'WPFE_data.mat'],'forecast','measurement','total_cap'); 
TSwind.forecast = zeros(365,96);
TSwind.production = zeros(365:96);
TSwind.delta_T = time_period;
TSwind.n_samples = length(forecast);

for a = 1:365
    TSwind.forecast(a,:) = forecast((a-1)*96+1:a*96);
    TSwind.production(a,:) = measurement((a-1)*96+1:a*96); 
end

% Load the demand data
load(['INPUT' filesep 'LOAD.mat'],'demand_net','other');

% the annual share of wind
share_wind = 0.30;

% rescale wind 
tot_wind = time_period*sum(sum(TSwind.production));
tot_load =  time_period*sum(sum(demand_net+other));
TSwind.cap_wind = share_wind/(tot_wind/tot_load);
TSwind.forecast = TSwind.forecast*TSwind.cap_wind;
TSwind.production = TSwind.production*TSwind.cap_wind;

% residual demand
res_demand = demand_net - TSwind.production;
cum_res_demand = sum(res_demand,2);

[min_res_demand,day_min_res_demand] = min(cum_res_demand);
[max_res_demand,day_max_res_demand] = max(cum_res_demand);

figure
plot(res_demand(day_min_res_demand,:))

figure
plot(res_demand(day_max_res_demand,:))

mean_res_demand = mean(res_demand);
diff_mean = res_demand - ones(365,1)*mean_res_demand;
cum_diff = sum(diff_mean.^2,2);

[useless,day_mean_dem] = min(cum_diff);

figure
plot(res_demand(day_mean_dem,:))

shift_res_demand = [zeros(365,1) res_demand];
shift_res_demand(1,1) = shift_res_demand(1,2);
for a = 2:365
    shift_res_demand(a,1) = shift_res_demand(a-1,end);
end
shift_res_demand(:,end) =[];
var_demand = res_demand - shift_res_demand;
cum_var_demand = sum(var_demand.^2,2);

cum_var_demand(32,:) = 0; 
cum_var_demand(109,:) = 0;
[useless,day_max_var] = max(cum_var_demand);

figure
plot(res_demand(day_max_var,:))
