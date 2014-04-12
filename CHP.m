clear all 
close all
clc
path(path,'C:\GAMS\win64\24.2')
clear gamso;
 
%% General information
% 
% Regarding the units
% * Prices in euros
% * Power in MW
% * Energy in MWh
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('%%%%%%%%%%%%%%%%%%%      CHP     %%%%%%%%%%%%%%%%%%%%%%%')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('          ')
 
%{
------------------------------------------------------------------
P_g= price electricity grid (spotprice)
P_c= price electricity consumption
P_i= price electricity imbalance
P_n= price natural gas
P_st= price startup

E_i0= electricity demand imbalance
E_b= electricity demand bidding

Q_H= heat demand house

Nb= thermal efficiency boiler
Ns= thermal efficiency storage
Ae= electrical efficiency CHP
Aq= thermal efficiency CHP

Pi_s= probability scenario

Ecap_lo= minimal electric energy supply CHP
Ecap_up= maximal electric energy supply CHP

Qcap_up= maximal thermal energy supply boiler

Cs= storage tank capacity
Cu= maximum number of CHP's
S= number of scenario's

x1= Heat output from the CHP 
X2= Auxiliary boiler output 
qst=Level of the heat storage tank 
qstout= Output of the heat storage 
%}


%---------------------------------------------
%% VARIABLES
%-----------------------------------
%LOADING VARIABLES 
%Input variables
%-----------------------------------

disp('Loading variables...')

load Spot2013 % spotprice of the electricity
load Imbalance2013 % imbalance prices
load Profiles2012 % heat and electricity demand of different houses, hotels and offices
%load Scenarios % different stochastic scenarios for wind imbalance
load WindTurbine2012 % Normalised wind turbine data (forecast and production)

%-----------------------------------
%DEFINING VARIABLES 
%Defining variables
%-----------------------------------

disp('Defining variables...')

% Input constants
%----------------

Ae0 = 0.4; % CHP electrical efficiency
Aq0 = 0.45; % CHP thermal efficiency
Nb0 = 0.98; % Boiler thermal efficiency
Ns0 = 0.998; % Storage heat efficiency

Cu_sh = 2; % shown chp
S_sh = 3; % shown scenario

% Input time (i)
%---------------

week = 3; % week of the year
sample_d = 2; % [days] sample duration

% Input units (n)
%----------------

N = 6; % number of units (CHP's, houses, etc.)

% Input scenarios (s)
%--------------------

S_gen = 100; % number of generated scenarios
S_red = 10; % number of reduced scenarios
S = S_red; % number of scenarios

%-----------------------------------
%INITIALISING VARIABLES 
%Input variables
%-----------------------------------

disp('Initialising variables...')

% Time (i)
%---------

% General
quarters = 4; % quarters per hour
hours = 24; % hours per day
days_w = 7; % days per week
days_y = 365; % days per year

timestep = 1/quarters; % [hours] 15 minutes

sample_h = sample_d*24; % [hours] sample duration
sample_q = sample_h * quarters; % [quarters] sample duration

period_d = week * days_w; % [days] period of the year
period_h = period_d * hours; % [hours] period of the year
period_q = period_h * quarters; % [quarters] period of the year

period_s = floor(period_d/sample_d); % number of complete sample durations in the previous period
year_s = floor(days_y/sample_d); % number of complete sample durations in a year
sample_pb = sample_q-(period_d/sample_d-period_s)+1*sample_q:sample_q; % uncomplete part of the sample (beginning)
sample_pe = 1:(days_y-year_s*sample_d)*hours*quarters-numel(sample_pb); % uncomplete part of the sample (end)

% Postprocessing
time_i = (1:sample_q)'; % time index vector sample
time_q = time_i/quarters; % [quarters] time vector sample

sample_i = [1:sample_q] + period_q; % sample period index

total_q=period_q+(sample_q-1); % [quarters] period before and during sample period

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

% Dimensioning of the CHP's (figure)
figure(1)
hold on
p1=plot(duration(:,Cu_sh),thermalload(:,Cu_sh));
h2=rectangle('Position',[0,0,duration_opt(Cu_sh),thermalload_opt(Cu_sh)]);
p2=plot(nan,nan,'s','markeredgecolor',get(h2,'edgecolor'));
legend([p1,p2],'Thermal load','Largest rectangle');
title(['Heat-load duration diagram for unit ',num2str(Cu_sh)]);
xlabel('Load duration [h]');
ylabel('Thermal load [MW_t]');

% Scenarios (s)
%--------------

% Actual production
actuals_y = TSwind.production; % Production for a 1MW wind turbine
actuals_s = actuals_y(sample_i); % Production for a 1MW wind turbine during the sample time

% Imbalance prices (relation to error)
price_actual_y = zeros(size(actuals_y));
tempPOS = POS;
tempNEG = -NEG;

price_actual_y(actuals_y>=0) = tempPOS(actuals_y>=0);
price_actual_y(actuals_y<0) = tempNEG(actuals_y<0);

price_actual_s = price_actual_y(sample_i,:);

price_posActual_y = POS; % [�/MWh] imbalance electricity price during the year
price_posActual_s = price_posActual_y(sample_i); % [�/MWh] imbalance electricity price during the sample duration

price_negActual_y = NEG; % [�/MWh] imbalance electricity price during the year
price_negActual_s = price_negActual_y(sample_i); % [�/MWh] imbalance electricity price during the sample duration

% Price
%------

% Gas
gasPrice = 0.040/2*1000; % [�/MWh] gas price (http://epp.eurostat.ec.europa.eu/statistics_explained/index.php/Electricity_and_natural_gas_price_statistics;http://appsso.eurostat.ec.europa.eu/nui/show.do?dataset=nrg_pc_204&lang=en;http://appsso.eurostat.ec.europa.eu/nui/show.do?dataset=nrg_pc_205&lang=en;http://appsso.eurostat.ec.europa.eu/nui/show.do?dataset=nrg_pc_202&lang=en;http://appsso.eurostat.ec.europa.eu/nui/show.do?dataset=nrg_pc_203&lang=en)
price_gas_y = gasPrice*ones(days_y*hours*quarters,1); % [�/MWh] gas price during the year
price_gas_s = price_gas_y(sample_i); % [�/MWh] gas price during the sample duration

% Consumption
price_elecC_d = reshape(repmat([0.15 0.15 0.15 0.15 0.15 0.15 0.15 0.22 0.22 0.22 0.22 0.22 0.22 0.22 0.22 0.22 0.22 0.22 0.22 0.22 0.22 0.22 0.15 0.15],quarters,1),hours*quarters,1)*1000; % [�/MWh] electricity price (day/night, one day)
price_elecC_y = repmat(price_elecC_d, days_y, 1); % [�/MWh] consumption electricity price during the year
price_elecC_s = price_elecC_y(sample_i); % [�/MWh] consumption electricity price during the sample duration

% Grid (spotprice)
price_elecS_y = reshape(repmat(spot,quarters,1),days_y*hours*quarters,1); % [�/MWh] spot electricity price during the year
price_elecS_s = price_elecS_y(sample_i); % [�/MWh] spot electricity price during the sample duration

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
        
elecD_s = elecD_y(sample_i,:); % [MWh] electric demand during the sample duration

% energy_y = reshape(repmat(mhouse1e+mhouse2e+mhouse3e,quarters,1),days_y*hours*quarters,1);
% energy_s = energy_y(sample_i);

% Input electric wind generation
%-------------------------------

elecWF_y = TSwind.forecast; % [MWh] forecast of electricity production with wind turbine 
elecWF_s = elecWF_y(sample_i); % [MWh] forecast of electricity production with wind turbine during the sample duration

elecWP_y = TSwind.production; % [MWh] electricity production with wind turbine 
elecWP_s = elecWP_y(sample_i); % [MWh] electricity production with wind turbine during the sample duration

erW_y = elecWP_y - elecWF_y; % [MWh] error between production and forecast
erW_s = erW_y(sample_i); % [MWh] error between production and forecast during the sample duration

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
        
heatD_s = heatD_y(sample_i,:); % [MWh] heat demand during the sample duration

% Capacity restrictions
% ---------------------

tank_cap = 2*max(heatD_y,[],1); % [MWh] Storage tank capacity

Ecap_upVar = 1/quarters*thermalload_opt/Aq0*Ae0;
Ecap_loVar = 0.5*Ecap_upVar;
Qcap_upVar = max(heatD_y,[],1)';

%-----------------------------------
%GENERATE FORECAST SCENARIOS
%generate the scenarios for the optimisation
%-----------------------------------

disp('Generating forecast scenarios...')

% Wind generation
%----------------

disp('  * Wind generation')

MU = zeros(1,numel(sample_i)); % average for random matrix
eps_init = 75; % initial value for covariance matrix
coef_fit_var = [-0.000028,0.0013,0.0037]; % coefficients of quadratic fit that compensates for the higher variability of wind at intermediate wind power forecasts 
norm_var = 0.0188; % normalization of quadratic fit
load('WPFE_ep.mat','pdf_stable','w_int','error_vec');
cdf_stable = cumsum(pdf_stable')';
load('WPFE_ee.mat','pdf_stable_ee','w_ee');
ff_mode = 'Old'; % scenario reduction mode (Old or Advanced)

[scen_lim,scen_red] = SCEN(S_gen,MU,eps_init,coef_fit_var,norm_var,elecWF_s,...
    sample_q,cdf_stable,pdf_stable,pdf_stable_ee,w_int,w_ee,error_vec,...
    S_red-1,erW_s,ff_mode);

% Spotprice
%----------

disp('  * Spotprice')

% TODO
% Based upon historic data?

% Heat demand
%------------

disp('  * Heat demand')

% TODO
% Based upon historic data?

% Imbalance prices
%-----------------

disp('  * Imbalance prices')

% TODO
% Based upon historic data?

% Gasprice
%---------

% TODO
% Probably not necessary since one price is used for the whole year

%-----------------------------------
%POSTPROCESSING SCENARIOS
%Postprocessing
%-----------------------------------

disp('Postprocessing scenarios...')

% Wind generation
%----------------

disp('  * Wind generation')

% Imbalance forecasts (error)
imbal_y = [scen_red.error(sample_pb); repmat(scen_red.error,year_s,1); scen_red.error(sample_pb)]; % Error is for wind turbine of 1 MW
imbal_s = scen_red.error;
Pi_st = scen_red.prob';

% Imbalance prices (relation to error)
price_imbal_y = zeros(size(imbal_y));
tempPOS = repmat(POS,1,S);
tempNEG = -repmat(NEG,1,S);

price_imbal_y(imbal_y>=0) = tempPOS(imbal_y>=0);
price_imbal_y(imbal_y<0) = tempNEG(imbal_y<0);

price_imbal_s = price_imbal_y(sample_i,:);

price_posImbal_y = POS; % [�/MWh] imbalance electricity price during the year
price_posImbal_s = price_posImbal_y(sample_i); % [�/MWh] imbalance electricity price during the sample duration

price_negImbal_y = NEG; % [�/MWh] imbalance electricity price during the year
price_negImbal_s = price_negImbal_y(sample_i); % [�/MWh] imbalance electricity price during the sample duration

% Spotprice
%----------

disp('  * Spotprice')

% TODO

% Heat demand
%------------

disp('  * Heat demand')

% TODO

% Imbalance prices
%-----------------

disp('  * Imbalance prices')

% TODO

%-----------------------------------
%PREPARE FOR DAY AHEAD OPTIMALISATION
%Input variables in structure form 
%-----------------------------------

disp('Preparing for day-ahead optimisation...')

[i,n,s,P_g,P_c,P_i,P_n,P_st,E_i0,E_i1,Q_H,Nb,Ns,Ae,Aq,Pi_s,...
    Ecap_lo,Ecap_up,Qcap_up,Cs,bid,bid_bool,dt] = GAMSWRITE(sample_q,N,...
    S,price_elecS_s,price_elecC_s,price_imbal_s,price_gas_s,imbal_s,...
    heatD_s,Nb0,Ns0,Ae0,Aq0,Pi_st,Ecap_loVar,Ecap_upVar,Qcap_upVar,...
    tank_cap,0,0,quarters);

%% CALL GAMS 
%-----------------------------------
%OPTIMIZATION
%Calls GAMS routine 
%-----------------------------------

disp('Optimising day-ahead...')

% Calculate investment costs (optimal)
%-------------------------------------

o2 = 0;
down = 1;

% while true
%     BUYst.val=[ones(1,Cust) zeros(1,Cu-Cust)];
%     wgdx('inputs', i,n,s,P_g,P_c,P_i,P_n,P_st,E_i0,E_i1,Q_H,Nb,Ns,Ae,Aq,Pi_s,Ecap_lo,Ecap_up,Qcap_up,Cs,dt,BUYst);
%     path(path,'C:\GAMS\win64\24.2')
%     gams('CHP');%LOCAL PRICE IS TAKEN INTO ACCOUNT 
%     
%     rs.name = 'obj';
%     r = rgdx ('results', rs);
%     obj=r.val(:,1);
%     o1 = obj/Cust;
%     [obj Cust o1]
%     
%     o1
%     break;
%     
%     if o1 >= o2
%         o2 = o1;
%         if down == 1
%             Cust = Cust - 1;
%             if Cust == 0; break; end;
%         else
%             Cust = Cust + 1;
%             if Cust > Cu; break; end;
%         end
%     else
%         if down == 1
%             Cust = Cust + 2;
%             if Cust > Cu; break; end;
%             down = 0;
%         else
%             break;
%         end
%     end
% end

% Custopt = Cust;

% Calculate investment costs 
%---------------------------

% Optimal is of course one CHP. Every CHP you add will have a slightly
% higher investment cost (imbalance remaining is less and thus more
% expensive to compensate). This part of the code let you calculate a
% random case (you can select one from the output of the previous part).

Cust = N; % Custopt;

BUYst.val=[ones(1,Cust) zeros(1,N-Cust)];
    wgdx('inputs', i,n,s,P_g,P_c,P_i,P_n,P_st,E_i0,E_i1,Q_H,Nb,Ns,Ae,Aq,Pi_s,Ecap_lo,Ecap_up,Qcap_up,Cs,bid,bid_bool,dt);

    gams('CHP'); % Day-ahead optimisation
    
    rs.name = 'obj';
    r = rgdx ('results', rs);
    obj=r.val(:,1);
    o1 = obj/Cust;

%% RESULTS 
%-----------------------------------
%LOADING RESULTS
%Reading parameters
%-----------------------------------

disp('Loading results day-ahead optimisation...')

[da.obj,da.R_b,da.R_ir,da.FC_bc,da.m_fCHP,da.m_fB,da.Q_CHP,da.Q_B,da.DeltaQ_S,da.Q_S,da.E_CHP,da.E_i,da.E_b,da.ON] = GAMSREAD(time_i,sample_q,N,S,sample_i);

%-----------------------------------
%POSTPROCESSING RESULTS
%Postproces data
%-----------------------------------

disp('Postprocessing results day-ahead optimisation...')

%R_b = E_b.*P_g.val;
%R_i = abs(E_i*Pi_st).*P_i.val;

profit = da.R_b + da.R_ir - da.FC_bc;

profit_t = sum(profit,1); % [�] total profit

%% Comparison
%-----------------------------------
%CALCULATE ACTUALS
%Calculate
%-----------------------------------

disp('Calculating actuals...')

%Input sets
%----------

i.uels = {{1:sample_q}};             
n.uels = {{1:N}};             
s.uels = {{1}};             


%Input parameters
%----------------

[i,n,s,P_g,P_c,P_i,P_n,P_st,E_i0,E_i1,Q_H,Nb,Ns,Ae,Aq,Pi_s,...
    Ecap_lo,Ecap_up,Qcap_up,Cs,bid,bid_bool,dt] = GAMSWRITE(sample_q,N,...
    S,price_elecS_s,price_elecC_s,price_actual_s,price_gas_s,actuals_s,...
    heatD_s,Nb0,Ns0,Ae0,Aq0,Pi_st,Ecap_loVar,Ecap_upVar,Qcap_upVar,...
    tank_cap,da.E_b(1),1,quarters);

BUYst.val=[ones(1,Cust) zeros(1,N-Cust)];
    wgdx('inputs', i,n,s,P_g,P_c,P_i,P_n,P_st,E_i0,E_i1,Q_H,Nb,Ns,Ae,Aq,Pi_s,Ecap_lo,Ecap_up,Qcap_up,Cs,bid,bid_bool,dt);

    gams('CHP'); % Calculate actuals

[a.obj,a.R_b,a.R_ir,a.FC_bc,a.m_fCHP,a.m_fB,a.Q_CHP,a.Q_B,a.DeltaQ_S,a.Q_S,a.E_CHP,a.E_i,a.E_b,a.ON] = GAMSREAD(time_i,sample_q,N,S,sample_i);

%-----------------------------------
%COMPARE WITH ACTUALS
%Comparison
%-----------------------------------

disp('Comparing actuals with day-ahead optimisation...')

% TODO

%% Figures
%-----------------------------------
%SHOWING RESULTS
%Making figures
%-----------------------------------

disp('Drawing figures...')

figure(2)
plot(time_q,[imbal_s(:,S_sh), da.E_b, da.E_i(:,S_sh)])
title(['Imbalance for scenario ',num2str(S_sh)]);
legend('Wind imbalance','Bidding','Imbalance reduction');

figure(3)
subplot(2,1,1)
plot(time_q,da.Q_S(:,:,S_sh))
title(['Storage for scenario ',num2str(S_sh), ' (all units)']);
subplot(2,1,2)
plot(time_q,da.DeltaQ_S(:,:,S_sh))
title(['Storage heat supply for scenario ',num2str(S_sh), ' (all units)']);

figure(4)
plot(time_q,[da.Q_S(:,Cu_sh,S_sh),da.DeltaQ_S(:,Cu_sh,S_sh)]);
legend('Storage','Storage heat supply');
title(['Storage unit ',num2str(Cu_sh),' for scenario ',num2str(S_sh)]);
xlabel('Time [h]');
ylabel('Heat [kWh]');

figure(5)
subplot(2,1,1)
plot(time_q,da.E_CHP(:,:,S_sh))
title(['CHP energy supply for scenario ',num2str(S_sh)]);
xlabel('Time [h]');
ylabel('Energy [kWh]');
subplot(2,1,2)
plot(time_q,da.Q_CHP(:,:,S_sh))
title(['CHP heat supply for scenario ',num2str(S_sh)]);
xlabel('Time [h]');
ylabel('Heat [kWh]');

figure(6)
plot(time_q,[heatD_s(:,Cu_sh),da.Q_B(:,Cu_sh,S_sh),da.Q_CHP(:,Cu_sh,S_sh)])
legend('House heat demand','Boiler heat supply','CHP heat supply');
title(['Heating unit ',num2str(Cu_sh),' for scenario ',num2str(S_sh)]);
xlabel('Time [h]');
ylabel('Heat [kWh]');

figure(7)
plot(time_q,[sum(da.m_fB(:,:,S_sh),2),sum(da.m_fCHP(:,:,S_sh),2)])
legend('Total fuelflow Boilers','Total fuelflow CHPs');
title(['Total fuelflow for scenario ',num2str(S_sh)]);
xlabel('Time [h]');
ylabel('Fuelflow [kWh]');

figure(8)
plot(time_q,[imbal_s(:,S_sh),da.E_CHP(:,Cu_sh,S_sh),da.E_b(:),da.E_i(:,S_sh)])
legend('Imbalance','CHP energy supply unit 1','Bidding','Total imbalance reduction');
title(['Energy unit ',num2str(Cu_sh),' for scenario ',num2str(S_sh)]);
xlabel('Time [h]');
ylabel('Energy [kWh]');

figure(9)
plot(time_q,[imbal_s*Pi_st,permute(sum(da.E_CHP,2),[1 3 2])*Pi_st,da.E_b(:),da.E_i*Pi_st])
legend('Imbalance','CHP energy supply','Bidding','Imbalance reduction');
title(['Energy for scenario ',num2str(S_sh)]);
xlabel('Time [h]');
ylabel('Energy [kWh]');

figure(10)
plot(time_q,[da.FC_bc,da.R_b,da.R_ir,profit]);
legend('Fuel consumption cost','Bidding revenue','Imbalance reduction revenue','profit');
title('Total profit');
xlabel('Time [h]');
ylabel('Money [�]');

% figure(1)
% subplot(2,1,1);
% plot(i,d1,'c')
% hold on 
% plot(i,CHP_thermal1,'r')
% hold on 
% plot(i,boiler1,'g')
% hold on
% plot(i,charge1)
% xlabel('Time [h]')
% legend('House heat demand', 'CHP heat supply', 'Boiler heat supply', 'Storage heat supply')

% subplot(2,1,2);
% plot(i,-charge1)
% hold on 
% plot(i,qst1,'k--')
% xlabel('Time [h]')
% legend('Heatflow to storage', 'Storage')

%% Data
%-----------------------------------
%SHOWING RESULTS
%show data
%-----------------------------------

disp('    ')
disp('DATA')
disp([' Profit: �',num2str(profit_t)])