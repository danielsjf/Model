function [Gen, Year] = INITYEAR(Par)
%-----------------------------------
%INITIALISING VARIABLES 
%Input variables
%-----------------------------------

AE0 = Par.AE0;                     
AQ0 = Par.AQ0;                     
Nb0 = Par.Nb0;                     
Ns0 = Par.Ns0;
MinCap = Par.MinCap;
E_bSingle = Par.E_bSingle;         
Yearstart = Par.Yearstart;         
week = Par.week;                   
sample_d = Par.sample_d;           
Dnom = Par.Dnom;                   
N = Par.N;                         
S_gen = Par.S_gen;                 
S_red = Par.S_red;

load Spot2013 % spotprice of the electricity
load Imbalance2013 % imbalance prices
load Profiles2012 % heat and electricity demand of different houses, hotels and offices
%load Scenarios % different stochastic scenarios for wind imbalance
load WindTurbine2012 % Normalised wind turbine data (forecast and production)
load CHPsets

disp('Initialising variables...')

% Time (i)
%---------

% General
quarters = 4; % quarters per hour
hours = 24; % hours per day
days_w = 7; % days per week
days_y = 365; % days per year

timestep = 1/quarters; % [hours] 15 minutes

sample_h = 2*24+(24-Par.Dnom); % [hours] sample duration
sample_q = sample_h * quarters; % [quarters] sample duration

% Postprocessing
time_i = (1:sample_q)'; % time index vector sample
time_q = time_i/quarters; % [quarters] time vector sample

% Units (n)
%----------

% General
N_data = 6; % number of units for which there is unique data
N_mult = ceil(N/N_data); % smallest multiple of N_data, larger than N (number of units)

% Dimensioning of the CHP's (pre-calculations)
bins = 200;
amount = zeros(N_data,bins);
thermalload = zeros(N_data,bins);

[amount(1,:), thermalload(1,:)] = hist(mhouse1h,bins); % Make histograms
[amount(2,:), thermalload(2,:)] = hist(mhouse2h,bins);
[amount(3,:), thermalload(3,:)] = hist(mhouse3h,bins);
[amount(4,:), thermalload(4,:)] = hist(mhouse4h,bins);
[amount(5,:), thermalload(5,:)] = hist(mhouse5h,bins);
[amount(6,:), thermalload(6,:)] = hist(mhouse6h,bins);

amount = repmat(amount,N_mult,1);
amount = amount(1:N,:);

thermalload = repmat(thermalload,N_mult,1);
thermalload = thermalload(1:N,:);

duration = cumsum(flipud(amount'),1);
thermalload = flipud(thermalload')/1000;

duration_opt = zeros(N,1); % [hours] Duration of the optimal thermal load in the heat-load diagram
thermalload_opt = zeros(N,1); % [MW] Optimal maximal thermal load of the CHP for the given load profile

% Dimensioning of the CHP's (real calculations)
for k=1:N
    loadsquares = duration(:,k).*thermalload(:,k); % Make an array of all the squares (duration x load)
    square_max = max(loadsquares); % Find the largest square
    j = find(loadsquares==square_max); % Find the ID of the largest square
    duration_opt(k) = duration(j,k); % Duration in the load diagram of the largest square
    thermalload_opt(k) = thermalload(j,k); % Load in the load diagram of the largest square
end

% Scenarios (s)
%--------------

% Actual production error
actuals_y = TSwind.production -  TSwind.forecast; % Production for a 1MW wind turbine

% Imbalance prices (relation to error)
price_actual_y = zeros(size(actuals_y));
tempPOS = POS;
tempNEG = -NEG;

price_actual_y(actuals_y>=0) = tempPOS(actuals_y>=0);
price_actual_y(actuals_y<0) = tempNEG(actuals_y<0);

price_posActual_y = POS; % [€/MWh] imbalance electricity price during the year

price_negActual_y = NEG; % [€/MWh] imbalance electricity price during the year

Pi_actual = 1; % The probability of the actual scenario

% Price
%------

% Gas
gasPrice = 0.040*1000; % [€/MWh] gas price (http://epp.eurostat.ec.europa.eu/statistics_explained/index.php/Electricity_and_natural_gas_price_statistics;http://appsso.eurostat.ec.europa.eu/nui/show.do?dataset=nrg_pc_204&lang=en;http://appsso.eurostat.ec.europa.eu/nui/show.do?dataset=nrg_pc_205&lang=en;http://appsso.eurostat.ec.europa.eu/nui/show.do?dataset=nrg_pc_202&lang=en;http://appsso.eurostat.ec.europa.eu/nui/show.do?dataset=nrg_pc_203&lang=en)
price_gas_y = gasPrice*ones(days_y*hours*quarters,1); % [€/MWh] gas price during the year

% Consumption
price_elecC_d = reshape(repmat([0.15 0.15 0.15 0.15 0.15 0.15 0.15 0.22 0.22 0.22 0.22 0.22 0.22 0.22 0.22 0.22 0.22 0.22 0.22 0.22 0.22 0.22 0.15 0.15],quarters,1),hours*quarters,1)*1000; % [€/MWh] electricity price (day/night, one day)
price_elecC_y = repmat(price_elecC_d, days_y, 1); % [€/MWh] consumption electricity price during the year

% Grid (spotprice)
price_elecS_y = reshape(repmat(spot,quarters,1),days_y*hours*quarters,1); % [€/MWh] spot electricity price during the year

% Input electric house demand
%----------------------------

elecD_y = [ reshape(repmat(mhouse1e',quarters,1),days_y*hours*quarters,1),...
            reshape(repmat(mhouse2e',quarters,1),days_y*hours*quarters,1),...
            reshape(repmat(mhouse3e',quarters,1),days_y*hours*quarters,1),...
            reshape(repmat(mhouse4e',quarters,1),days_y*hours*quarters,1),...
            reshape(repmat(mhouse5e',quarters,1),days_y*hours*quarters,1),...
            reshape(repmat(mhouse6e',quarters,1),days_y*hours*quarters,1)]/1000; % [MWh] electric demand during the year
        
elecD_y = repmat(elecD_y,1,N_mult);
elecD_y = elecD_y(:,1:N);

% energy_y = reshape(repmat(mhouse1e+mhouse2e+mhouse3e,quarters,1),days_y*hours*quarters,1);
% energy_s = energy_y(sample_i);

% Input electric wind generation
%-------------------------------

elecWF_y = TSwind.forecast; % [MWh] forecast of electricity production with wind turbine 

elecWP_y = TSwind.production; % [MWh] electricity production with wind turbine 

erW_y = elecWP_y - elecWF_y; % [MWh] error between production and forecast

% Input heat demand
%------------------

heatD_y = [ reshape(repmat(mhouse1h',quarters,1),days_y*hours*quarters,1),...
            reshape(repmat(mhouse2h',quarters,1),days_y*hours*quarters,1),...
            reshape(repmat(mhouse3h',quarters,1),days_y*hours*quarters,1),...
            reshape(repmat(mhouse4h',quarters,1),days_y*hours*quarters,1),...
            reshape(repmat(mhouse5h',quarters,1),days_y*hours*quarters,1),...
            reshape(repmat(mhouse6h',quarters,1),days_y*hours*quarters,1)]/1000; % [MWh] heat demand during the year
        
heatD_y = repmat(heatD_y,1,N_mult);
heatD_y = heatD_y(:,1:N);

% Capacity restrictions
% ---------------------

tank_cap = 2*max(heatD_y,[],1); % [MWh] Storage tank capacity

CHPData = Set{Par.CHPset};

Ecap_upVar = 1/quarters*[CHPData(3).EC, CHPData(1).EC, CHPData(2).EC, CHPData(1).EC, CHPData(1).EC, CHPData(2).EC]/1000;
Ecap_loVar = MinCap*Ecap_upVar;
Qcap_upVar = max(heatD_y,[],1)';

%% Gen

Gen.quarters = quarters;        
Gen.hours = hours;              
Gen.days_w = days_w;            
Gen.days_y = days_y;            
Gen.timestep = timestep;        
Gen.sample_h = sample_h;        
Gen.sample_q = sample_q;        
Gen.time_i = time_i;            
Gen.time_q = time_q;            
Gen.tank_cap = tank_cap;        
Gen.Ecap_upVar = Ecap_upVar;    
Gen.Ecap_loVar = Ecap_loVar;    
Gen.Qcap_upVar = Qcap_upVar;

%% Year

Year.thermalload = thermalload;                
Year.duration_opt = duration_opt;              
Year.thermalload_opt = thermalload_opt;        
Year.actuals_y = actuals_y;                    
Year.price_actual_y = price_actual_y;          
Year.tempPOS = tempPOS;                        
Year.tempNEG = tempNEG;                        
Year.price_posActual_y = price_posActual_y;    
Year.price_negActual_y = price_negActual_y;    
Year.Pi_actual = Pi_actual;                                           
Year.gasPrice = gasPrice;                      
Year.price_gas_y = price_gas_y;                
Year.price_elecC_y = price_elecC_y;            
Year.price_elecS_y = price_elecS_y;            
Year.elecD_y = elecD_y;                        
Year.elecWF_y = elecWF_y;                      
Year.elecWP_y = elecWP_y;                      
Year.erW_y = erW_y;                            
Year.heatD_y = heatD_y; 
end