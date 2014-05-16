clear all 
close all
clc
path(path,'C:\GAMS\win64\24.2')
path(path,[pwd filesep 'DataHash'])
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

E_i= electricity demand imbalance
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

Par.spot = DataHash(spot); % Hash of the spotprice data to check whether the data has been changed
Par.imbalance = DataHash([MDP, MIP, NEG, NRV, POS, SI, alpha]); % Hash of the imbalance data to check whether the data has been changed
Par.profiles = DataHash([mhouse1h, mhouse2h, mhouse3h, mhouse4h, mhouse5h, mhouse6h]); % Hash of the heat profile data to check whether the data has been changed
Par.windturbine = DataHash(TSwind); % Hash of the windturbine data to check whether the data has been changed


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

RunAgain = 1; % 0 is no run as the parameters already were used in the past, 1 is a run as the parameters are already used

E_bSingle = 0; % 1 is single bid for the whole period, 0 is one bid per quarter

% Input time (i)
%---------------

Yearstart = 2; % First day of 2013 was a Tuesday
week = [3, 16, 28]; % week of the year
sample_d = 7; % [days] sample duration (every time a full week)
Dnom = 12; % [hours] hour of the nomination deadline

% Input units (n)
%----------------

N = 6; % number of units (CHP's, houses, etc.)

% Input scenarios (s)
%--------------------

S_gen = 5; % number of generated scenarios
S_red = 2; % number of reduced scenarios
S = S_red; % number of scenarios

% Convert
%--------

Par.Ae0 = Ae0;                    
Par.Aq0 = Aq0;                    
Par.Nb0 = Nb0;                    
Par.Ns0 = Ns0;                    
Par.E_bSingle = E_bSingle;        
Par.Yearstart = Yearstart;        
Par.week = week;                  
Par.sample_d = sample_d;          
Par.Dnom = Dnom;                  
Par.N = N;                        
Par.S_gen = S_gen;                
Par.S_red = S_red; 

% Hash
%-----

parHash = DataHash(Par); % Hash of the parameters data to check whether the data has been changed
if exist('Results.mat', 'file') == 2
    load('Results.mat');
    parHashOld = cell(size(resultsOld,2),1);
    for k = 1:size(resultsOld,2)
        parHashOld(k) = cellstr(resultsOld(k).parHash);
    end
else
    parHashOld = '';
end

nR = find(strcmp(parHash,parHashOld)); % Find if there is already an other version with the same calculations

Calc = [];

if size(nR,1) == 1
    if RunAgain == 0
        Calc = 0;
    else
        Calc = 1;
    end
else
    Calc = 1;
end

% %-----------------------------------
% %INITIALISING VARIABLES 
% %Input variables
% %-----------------------------------

[Gen, Year] = INITYEAR(Par);

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
CHPBool = Year.CHPBool;                         
gasPrice = Year.gasPrice;                       
price_gas_y = Year.price_gas_y;                 
price_elecC_y = Year.price_elecC_y;             
price_elecS_y = Year.price_elecS_y;             
elecD_y = Year.elecD_y;                         
elecWF_y = Year.elecWF_y;                       
elecWP_y = Year.elecWP_y;                       
erW_y = Year.erW_y;                             
heatD_y = Year.heatD_y; 

aT.R_b = [];
aT.R_ir = [];
aT.FC_bc = [];
aT.m_fCHP = [];
aT.m_fB = [];
aT.Q_CHP = [];
aT.Q_B = [];
aT.DeltaQ_S = [];
aT.Q_S = [];
aT.E_CHP = [];
aT.E_ir = [];
aT.E_b = [];
aT.ON = [];

for xi = 1:numel(week) % sample week of the year
    for xj = 1:sample_d+1 % day of the sample week
%-----------------------------------
%INITIALISING VARIABLES 
%Input variables
%-----------------------------------

[Sample] = INITSAMPLE(Par,Gen,Year,xi,xj);

sample_i = Sample.sample_i;    
sample_ii = Sample.sample_ii;   
sample_pe = Sample.sample_pe;                     
total_q = Sample.total_q;                         
actuals_s = Sample.actuals_s;                     
price_actual_s = Sample.price_actual_s;           
price_posActual_s = Sample.price_posActual_s;     
price_negActual_s = Sample.price_negActual_s;     
price_gas_s = Sample.price_gas_s;                 
price_elecC_s = Sample.price_elecC_s;             
price_elecS_s = Sample.price_elecS_s;             
elecD_s = Sample.elecD_s;                         
elecWF_s = Sample.elecWF_s;                       
elecWP_s = Sample.elecWP_s;                       
erW_s = Sample.erW_s;                             
heatD_s = Sample.heatD_s;

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

if Calc == 1
    [scen_lim,scen_red] = SCEN(Par.S_gen,MU,eps_init,coef_fit_var,norm_var,elecWF_s,...
        sample_q,cdf_stable,pdf_stable,pdf_stable_ee,w_int,w_ee,error_vec,...
        Par.S_red-1,erW_s,ff_mode);
else
    Scen = resultsOld(nR).Scen;
    scen_lim = Scen.lim;
    scen_red = Scen.red;
end

Scen.lim = scen_lim;
Scen.red = scen_red;

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
imbal_y = zeros(size(heatD_y,1),S); %[scen_red.error(sample_pb); repmat(scen_red.error,year_s,1); scen_red.error(sample_pb)]; % Error is for wind turbine of 1 MW
imbal_s = scen_red.error;
Pi_st = scen_red.prob';

% Imbalance prices (relation to error)
price_imbal_y = zeros(size(imbal_y));
tempPOS = repmat(POS,1,S);
tempNEG = -repmat(NEG,1,S);

price_imbal_y(imbal_y>=0) = tempPOS(imbal_y>=0);
price_imbal_y(imbal_y<0) = tempNEG(imbal_y<0);

price_imbal_s = price_imbal_y(sample_i,:);

price_posImbal_y = POS; % [€/MWh] imbalance electricity price during the year
price_posImbal_s = price_posImbal_y(sample_i); % [€/MWh] imbalance electricity price during the sample duration

price_negImbal_y = NEG; % [€/MWh] imbalance electricity price during the year
price_negImbal_s = price_negImbal_y(sample_i); % [€/MWh] imbalance electricity price during the sample duration

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

[i,n,s,P_g,P_c,P_i,P_n,P_st,E_i,E_is,Q_H,Nb,Ns,Ae,Aq,Pi_s,...
    Ecap_lo,Ecap_up,Qcap_up,Cs,bid,bid_bool,bid_single,CHP_bool,dt] = GAMSWRITE(sample_q,N,...
    S,price_elecS_s,price_elecC_s,price_imbal_s,price_gas_s,imbal_s,...
    heatD_s,Par.Nb0,Par.Ns0,Par.Ae0,Par.Aq0,Pi_st,Ecap_loVar,Ecap_upVar,Qcap_upVar,...
    tank_cap,zeros(sample_q,1),0,Par.E_bSingle,CHPBool,quarters);

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
%     wgdx('inputs', i,n,s,P_g,P_c,P_i,P_n,P_st,E_i,E_is,Q_H,Nb,Ns,Ae,Aq,Pi_s,Ecap_lo,Ecap_up,Qcap_up,Cs,dt,BUYst);
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
    wgdx('inputs', i,n,s,P_g,P_c,P_i,P_n,P_st,E_i,E_is,Q_H,Nb,Ns,Ae,Aq,Pi_s,Ecap_lo,Ecap_up,Qcap_up,Cs,bid,bid_bool,bid_single,CHP_bool,dt);

if Calc == 1
    gams('CHP'); % Day-ahead optimisation
end

%% RESULTS 
%-----------------------------------
%LOADING RESULTS
%Reading parameters
%-----------------------------------

disp('Loading results day-ahead optimisation...')

if Calc == 1
    da.name = 'day-ahead optimisation';
    [da.obj,da.R_b,da.R_ir,da.FC_bc,da.m_fCHP,da.m_fB,da.Q_CHP,da.Q_B,da.DeltaQ_S,da.Q_S,da.E_CHP,da.E_ir,da.E_b,da.ON] = GAMSREAD(time_i,sample_q,N,S,sample_i);
else
    da = resultsOld(nR).da;
end

%-----------------------------------
%POSTPROCESSING RESULTS
%Postproces data
%-----------------------------------

disp('Postprocessing results day-ahead optimisation...')

%R_b = E_b.*P_g.val;
%R_i = abs(E_ir*Pi_st).*P_i.val;

profit = da.R_b + da.R_ir - da.FC_bc;

profit_t = sum(profit,1); % [€] total profit

%% Comparison
%-----------------------------------
%CALCULATE ACTUALS
%Calculate
%-----------------------------------

disp('Calculating actuals...')

% Input actuals
%--------------

S = 1;

% Input sets
%-----------

i.uels = {{1:sample_q}};             
n.uels = {{1:N}};             
s.uels{1,1} = 1:S;             


% Input parameters
%-----------------

% Actuals (with CHP; optimal)
[i,n,s,P_g,P_c,P_i,P_n,P_st,E_i,E_is,Q_H,Nb,Ns,Ae,Aq,Pi_s,...
    Ecap_lo,Ecap_up,Qcap_up,Cs,bid,bid_bool,bid_single,CHP_bool,dt] = GAMSWRITE(sample_q,N,...
    S,price_elecS_s,price_elecC_s,price_actual_s,price_gas_s,actuals_s,...
    heatD_s,Par.Nb0,Par.Ns0,Par.Ae0,Par.Aq0,Pi_actual,Ecap_loVar,Ecap_upVar,Qcap_upVar,...
    tank_cap,zeros(sample_q,1),0,Par.E_bSingle,CHPBool,quarters);

BUYst.val=[ones(1,Cust) zeros(1,N-Cust)];
    wgdx('inputs', i,n,s,P_g,P_c,P_i,P_n,P_st,E_i,E_is,Q_H,Nb,Ns,Ae,Aq,Pi_s,Ecap_lo,Ecap_up,Qcap_up,Cs,bid,bid_bool,bid_single,CHP_bool,dt);

if Calc == 1
    gams('CHP'); % Calculate actuals
end

if Calc == 1
    a_OPT.name = 'actuals (with CHP; optimal)';
    [a_OPT.obj,a_OPT.R_b,a_OPT.R_ir,a_OPT.FC_bc,a_OPT.m_fCHP,a_OPT.m_fB,a_OPT.Q_CHP,a_OPT.Q_B,a_OPT.DeltaQ_S,a_OPT.Q_S,a_OPT.E_CHP,a_OPT.E_ir,a_OPT.E_b,a_OPT.ON] = GAMSREAD(time_i,sample_q,N,S,sample_i);
else
    a_OPT = resultsOld(nR).a_OPT;
end

% Actuals (with CHP; day-ahead bidding)
[i,n,s,P_g,P_c,P_i,P_n,P_st,E_i,E_is,Q_H,Nb,Ns,Ae,Aq,Pi_s,...
    Ecap_lo,Ecap_up,Qcap_up,Cs,bid,bid_bool,bid_single,CHP_bool,dt] = GAMSWRITE(sample_q,N,...
    S,price_elecS_s,price_elecC_s,price_actual_s,price_gas_s,actuals_s,...
    heatD_s,Par.Nb0,Par.Ns0,Par.Ae0,Par.Aq0,Pi_actual,Ecap_loVar,Ecap_upVar,Qcap_upVar,...
    tank_cap,da.E_b,1,Par.E_bSingle,CHPBool,quarters);

BUYst.val=[ones(1,Cust) zeros(1,N-Cust)];
    wgdx('inputs', i,n,s,P_g,P_c,P_i,P_n,P_st,E_i,E_is,Q_H,Nb,Ns,Ae,Aq,Pi_s,Ecap_lo,Ecap_up,Qcap_up,Cs,bid,bid_bool,bid_single,CHP_bool,dt);

if Calc == 1
    gams('CHP'); % Calculate actuals
end    

if Calc == 1
    a.name = 'actuals (with CHP; day-ahead bidding)';
    [a.obj,a.R_b,a.R_ir,a.FC_bc,a.m_fCHP,a.m_fB,a.Q_CHP,a.Q_B,a.DeltaQ_S,a.Q_S,a.E_CHP,a.E_ir,a.E_b,a.ON] = GAMSREAD(time_i,sample_q,N,S,sample_i);
else
    a = resultsOld(nR).a;
end

% Actuals (without CHP; day-ahead bidding/optimal)
CHPBool = 0;
[i,n,s,P_g,P_c,P_i,P_n,P_st,E_i,E_is,Q_H,Nb,Ns,Ae,Aq,Pi_s,...
    Ecap_lo,Ecap_up,Qcap_up,Cs,bid,bid_bool,bid_single,CHP_bool,dt] = GAMSWRITE(sample_q,N,...
    S,price_elecS_s,price_elecC_s,price_actual_s,price_gas_s,actuals_s,...
    heatD_s,Par.Nb0,Par.Ns0,Par.Ae0,Par.Aq0,Pi_actual,Ecap_loVar,Ecap_upVar,Qcap_upVar,...
    tank_cap,zeros(sample_q,1),0,CHPBool,Par.E_bSingle,quarters);

BUYst.val=[ones(1,Cust) zeros(1,N-Cust)];
    wgdx('inputs', i,n,s,P_g,P_c,P_i,P_n,P_st,E_i,E_is,Q_H,Nb,Ns,Ae,Aq,Pi_s,Ecap_lo,Ecap_up,Qcap_up,Cs,bid,bid_bool,bid_single,CHP_bool,dt);

if Calc == 1
    gams('CHP'); % Calculate actuals
end

if Calc == 1
    a_NCHP.name = 'actuals (without CHP; day-ahead bidding/optimal)';
    [a_NCHP.obj,a_NCHP.R_b,a_NCHP.R_ir,a_NCHP.FC_bc,a_NCHP.m_fCHP,a_NCHP.m_fB,a_NCHP.Q_CHP,a_NCHP.Q_B,a_NCHP.DeltaQ_S,a_NCHP.Q_S,a_NCHP.E_CHP,a_NCHP.E_ir,a_NCHP.E_b,a_NCHP.ON] = GAMSREAD(time_i,sample_q,N,S,sample_i);
else
    a_NCHP = resultsOld(nR).a_NCHP;
end

%-----------------------------------
%SUM UP
%Sum up
%-----------------------------------

aT.R_b = [aT.R_b; a.R_b(sample_ii)];
aT.R_ir = [aT.R_ir; a.R_ir(sample_ii)];
aT.FC_bc = [aT.FC_bc; a.FC_bc(sample_ii)];
aT.m_fCHP = [aT.m_fCHP; a.m_fCHP(sample_ii,:)];
aT.m_fB = [aT.m_fB; a.m_fB(sample_ii,:)];
aT.Q_CHP = [aT.Q_CHP; a.Q_CHP(sample_ii,:)];
aT.Q_B = [aT.Q_B; a.Q_B(sample_ii,:)];
aT.DeltaQ_S = [aT.DeltaQ_S; a.DeltaQ_S(sample_ii,:)];
aT.Q_S = [aT.Q_S; a.Q_S(sample_ii,:)];
aT.E_CHP = [aT.E_CHP; a.E_CHP(sample_ii,:)];
aT.E_ir = [aT.E_ir; a.E_ir(sample_ii)];
aT.E_b = [aT.E_b; a.E_b(sample_ii)];
aT.ON = [aT.ON; a.ON(sample_ii,:)];

%-----------------------------------
%COMPARE WITH ACTUALS
%Comparison
%-----------------------------------

disp('Comparing actuals with day-ahead optimisation...')

% Actuals (with CHP; optimal)
%[a_OPT.GEN, a_OPT.CHP_data] = STATSSCEN(a_OPT,1,1);

% Actuals (with CHP; day-ahead bidding)
%[a.GEN, a.CHP_data] = STATSSCEN(a,1,1); 

% Actuals (without CHP; day-ahead bidding/optimal)
%[a_NCHP.GEN, a_NCHP.CHP_data] = STATSSCEN(a_NCHP,1,1);

    end
end

%% Figures
%-----------------------------------
%SHOWING RESULTS
%Making figures
%-----------------------------------

disp('Drawing figures...')

%DRAWFIG(da,time_q,imbal_s,heatD_s,Pi_st,Cu_sh,S_sh);

%% Data
%-----------------------------------
%SHOWING RESULTS
%show data
%-----------------------------------

disp('    ')
disp('DATA')

%% Export Data
%-----------------------------------
%Export data
%export data
%-----------------------------------

disp('    ')
disp('Export data')

results.Par = Par;
results.parHash = parHash;
results.da = da;
results.a = a;
results.a_OPT = a_OPT;
results.a_NCHP = a_NCHP;
results.Scen = Scen;
if Calc == 1 && exist('Results.mat', 'file') == 2
    resultsOld = [resultsOld results];
    save('Results.mat', 'resultsOld');
elseif Calc == 1
    resultsOld = results;
    save('Results.mat', 'resultsOld');
end