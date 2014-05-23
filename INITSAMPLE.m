function [Sample] = INITSAMPLE(Par,Gen,Year,xi,xj)
%-----------------------------------
%INITIALISING VARIABLES 
%Input variables
%-----------------------------------

%Par
spot = Par.spot;                   
imbalance = Par.imbalance;         
profiles = Par.profiles;           
windturbine = Par.windturbine;     
AE0 = Par.AE0;                     
AQ0 = Par.AQ0;                     
Nb0 = Par.Nb0;                     
Ns0 = Par.Ns0;                     
E_bSingle = Par.E_bSingle;         
Yearstart = Par.Yearstart;         
week = Par.week;                   
sample_d = Par.sample_d;           
Dnom = Par.Dnom;                   
N = Par.N;                         
S_gen = Par.S_gen;                 
S_red = Par.S_red;  

%Gen
quarters = Gen.quarters;         
hours = Gen.hours;               
days_w = Gen.days_w;             
days_y = Gen.days_y;             
timestep = Gen.timestep;         
sample_h = Gen.sample_h;         
sample_q = Gen.sample_q;         
time_i = Gen.time_i;             
time_q = Gen.time_q;             
tank_cap = Gen.tank_cap;         
Ecap_upVar = Gen.Ecap_upVar;     
Ecap_loVar = Gen.Ecap_loVar;     
Qcap_upVar = Gen.Qcap_upVar;

%Year
thermalload = Year.thermalload;                 
duration_opt = Year.duration_opt;               
thermalload_opt = Year.thermalload_opt;         
actuals_y = Year.actuals_y;                     
price_actual_y = Year.price_actual_y;           
tempPOS = Year.tempPOS;                         
tempNEG = Year.tempNEG;                         
price_posActual_y = Year.price_posActual_y;     
price_negActual_y = Year.price_negActual_y;     
Pi_actual = Year.Pi_actual;                                        
gasPrice = Year.gasPrice;                       
price_gas_y = Year.price_gas_y;                 
price_elecC_y = Year.price_elecC_y;             
price_elecS_y = Year.price_elecS_y;             
elecD_y = Year.elecD_y;                         
elecWF_y = Year.elecWF_y;                       
elecWP_y = Year.elecWP_y;                       
erW_y = Year.erW_y;                             
heatD_y = Year.heatD_y;  

% General

period_d = (week(xi)-1) * days_w + 6 - Par.Yearstart + xj; % [days] period of the year
period_h = period_d * hours + Dnom; % [hours] period of the year
period_q = period_h * quarters; % [quarters] period of the year

period_s = floor(period_d/Par.sample_d); % number of complete sample durations in the previous period
year_s = floor(days_y/Par.sample_d); % number of complete sample durations in a year
sample_pb = sample_q-(period_d/Par.sample_d-period_s)+1*sample_q:sample_q; % uncomplete part of the sample (beginning)
sample_pe = 1:(days_y-year_s*Par.sample_d)*hours*quarters-numel(sample_pb); % uncomplete part of the sample (end)

% Postprocessing

sample_i = [1:sample_q] + period_q; % sample period index

total_q=period_q+(sample_q-1); % [quarters] period before and during sample period

sample_ii = (24-Dnom)*quarters+1:sample_q-24*quarters;

% Units (n)
%----------

% Scenarios (s)
%--------------

% Actual production error
actuals_s = actuals_y(sample_i); % Production for a 1MW wind turbine during the sample time

% Imbalance prices (relation to error)
price_actual_s = price_actual_y(sample_i,:);

price_posActual_s = price_posActual_y(sample_i); % [€/MWh] imbalance electricity price during the sample duration

price_negActual_s = price_negActual_y(sample_i); % [€/MWh] imbalance electricity price during the sample duration

% Price
%------

% Gas
price_gas_s = price_gas_y(sample_i); % [€/MWh] gas price during the sample duration

% Consumption
price_elecC_s = price_elecC_y(sample_i); % [€/MWh] consumption electricity price during the sample duration

% Grid (spotprice)
price_elecS_s = price_elecS_y(sample_i); % [€/MWh] spot electricity price during the sample duration

% Input electric house demand
%----------------------------
        
elecD_s = elecD_y(sample_i,:); % [MWh] electric demand during the sample duration

% energy_y = reshape(repmat(mhouse1e+mhouse2e+mhouse3e,quarters,1),days_y*hours*quarters,1);
% energy_s = energy_y(sample_i);

% Input electric wind generation
%-------------------------------

elecWF_s = elecWF_y(sample_i); % [MWh] forecast of electricity production with wind turbine during the sample duration

elecWP_s = elecWP_y(sample_i); % [MWh] electricity production with wind turbine during the sample duration

erW_s = erW_y(sample_i); % [MWh] error between production and forecast during the sample duration

% Input heat demand
%------------------
        
heatD_s = heatD_y(sample_i,:); % [MWh] heat demand during the sample duration

Sample.sample_i = sample_i;
Sample.sample_ii = sample_ii;
Sample.sample_pe = sample_pe;
Sample.total_q = total_q;
Sample.actuals_s = actuals_s;                    
Sample.price_actual_s = price_actual_s;          
Sample.price_posActual_s = price_posActual_s;    
Sample.price_negActual_s = price_negActual_s;    
Sample.price_gas_s = price_gas_s;                
Sample.price_elecC_s = price_elecC_s;            
Sample.price_elecS_s = price_elecS_s;            
Sample.elecD_s = elecD_s;                        
Sample.elecWF_s = elecWF_s;                      
Sample.elecWP_s = elecWP_s;                                                
Sample.heatD_s = heatD_s; 
end