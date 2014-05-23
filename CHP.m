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
load CHPsets % CHP data

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

CHPset = 3; % 1 = ICE; 2 = SOFC; 3 = PEM
CHPData = Set{CHPset};

AE0 = [CHPData(3).EE, CHPData(1).EE, CHPData(2).EE, CHPData(1).EE, CHPData(1).EE, CHPData(2).EE]/100; % CHP electrical efficiency
AQ0 = [CHPData(3).TE, CHPData(1).TE, CHPData(2).TE, CHPData(1).TE, CHPData(1).TE, CHPData(2).TE]/100; % CHP thermal efficiency
AE1 = [0 0 0 0 0 0]; % CHP electrical efficiency constant
AQ1 = [0 0 0 0 0 0]; % CHP thermal efficiency constant
Nb0 = 0.98; % Boiler thermal efficiency
Ns0 = 0.998; % Storage heat efficiency
MinCap = 0.3;

Cu_sh = 2; % shown chp
S_sh = 3; % shown scenario

RunAgain = 1; % 0 is no run as the parameters already were used in the past, 1 is a run as the parameters are already used

E_bSingle = 0; % 1 is single bid for the whole period, 0 is one bid per quarter

% Input time (i)
%---------------

Yearstart = 2; % First day of 2013 was a Tuesday
week = [2, 15, 28]; % week of the year
sample_d = 7; % [days] sample duration (every time a full week)
Dnom = 12; % [hours] hour of the nomination deadline

% Input units (n)
%----------------

N = 6; % number of units (CHP's, houses, etc.)

% Input scenarios (s)
%--------------------

S_gen = 100; % number of generated scenarios
S_red = 10; % number of reduced scenarios
S = S_red; % number of scenarios

% Convert
%--------

Par.CHPset = CHPset;
Par.AE0 = AE0;                    
Par.AQ0 = AQ0;   
Par.AE1 = AE1;                    
Par.AQ1 = AQ1;
Par.Nb0 = Nb0;                    
Par.Ns0 = Ns0;
Par.MinCap = MinCap;
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

Par2 = Par;
Par2.Dnom = 24;
[Gen2, ~] = INITYEAR(Par2);

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

SampleDAT = [];
SampleDAT = INITSTRUCT(SampleDAT,'Sample');

SampleDATT = [];
SampleDATT = INITSTRUCT(SampleDATT,'Sample');

SampleAT = [];
SampleAT = INITSTRUCT(SampleAT,'Sample');

SampleATT = [];
SampleATT = INITSTRUCT(SampleATT,'Sample');

aT = [];
aT = INITSTRUCT(aT,'CHPdata');

aTT = [];
aTT = INITSTRUCT(aTT,'CHPdata');

daT = [];
daT = INITSTRUCT(daT,'CHPdata');

daTT = [];
daTT = INITSTRUCT(daTT,'CHPdata');

for xi = 1:numel(week) % sample week of the year
    Q_Ss = zeros(N,1);
    on0 = zeros(N,1);
    UpPenal = 1; DoPenal = 10;
    if xi == 1; UpPenal = 1; DoPenal = 10; end
    for xj = 1:sample_d+1 % day of the sample week
%-----------------------------------
%INITIALISING VARIABLES 
%Input variables
%-----------------------------------

disp(['Day ', num2str(xj), ' of week ', num2str(xi)])

if xj ~= 1; Q_Ss = Q_SsDnom; end
if xj ~= 1; on0 = on0Dnom; end

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
heatD_s = Sample.heatD_s;

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

S = S_red;

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
        Par.S_red-1,actuals_s,ff_mode);
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

SampleDAT.sample_i = [SampleDAT.sample_i, sample_i(sample_ii)];
SampleDAT.imbal_s = [SampleDAT.imbal_s; imbal_s(sample_ii,:)];
SampleDAT.Pi_st = [SampleDAT.Pi_st; Pi_st'];
SampleDAT.actuals_s = [SampleDAT.actuals_s; actuals_s(sample_ii)];
SampleDAT.price_actual_s = [SampleDAT.price_actual_s; price_actual_s(sample_ii)];
SampleDAT.price_posActual_s = [SampleDAT.price_posActual_s; price_posActual_s(sample_ii)];
SampleDAT.price_negActual_s = [SampleDAT.price_negActual_s; price_negActual_s(sample_ii)];
SampleDAT.price_gas_s = [SampleDAT.price_gas_s; price_gas_s(sample_ii)];
SampleDAT.price_elecC_s = [SampleDAT.price_elecC_s; price_elecC_s(sample_ii)];
SampleDAT.price_elecS_s = [SampleDAT.price_elecS_s; price_elecS_s(sample_ii)];
SampleDAT.elecD_s = [SampleDAT.elecD_s; elecD_s(sample_ii,:)];
SampleDAT.elecWF_s = [SampleDAT.elecWF_s; elecWF_s(sample_ii)];
SampleDAT.elecWP_s = [SampleDAT.elecWP_s; elecWP_s(sample_ii)];
SampleDAT.heatD_s = [SampleDAT.heatD_s; heatD_s(sample_ii,:)];

SampleDATT.sample_i = [SampleDATT.sample_i, sample_i];
SampleDATT.imbal_s = [SampleDATT.imbal_s; imbal_s];
SampleDATT.Pi_st = [SampleDATT.Pi_st; Pi_st'];
SampleDATT.actuals_s = [SampleDATT.actuals_s; actuals_s];
SampleDATT.price_actual_s = [SampleDATT.price_actual_s; price_actual_s];
SampleDATT.price_posActual_s = [SampleDATT.price_posActual_s; price_posActual_s];
SampleDATT.price_negActual_s = [SampleDATT.price_negActual_s; price_negActual_s];
SampleDATT.price_gas_s = [SampleDATT.price_gas_s; price_gas_s];
SampleDATT.price_elecC_s = [SampleDATT.price_elecC_s; price_elecC_s];
SampleDATT.price_elecS_s = [SampleDATT.price_elecS_s; price_elecS_s];
SampleDATT.elecD_s = [SampleDATT.elecD_s; elecD_s];
SampleDATT.elecWF_s = [SampleDATT.elecWF_s; elecWF_s];
SampleDATT.elecWP_s = [SampleDATT.elecWP_s; elecWP_s];
SampleDATT.heatD_s = [SampleDATT.heatD_s; heatD_s]; 

[i,n,s,P_g,P_c,P_i,P_n,P_st,E_i,E_is,Q_H,Nb,Ns,Ae0,Aq0,Ae1,Aq1,Pi_s,...
    Ecap_lo,Ecap_up,Qcap_up,Cs,bid,bid_bool,bid_single,CHP_bool,QSs,ON0,UpPen,DoPen,Dur,dt] = GAMSWRITE(sample_q,N,...
    S,price_elecS_s,price_elecC_s,price_imbal_s,price_gas_s,imbal_s,...
    heatD_s,Nb0,Ns0,AE0,AQ0,AE1,AQ1,Pi_st,Ecap_loVar,Ecap_upVar,Qcap_upVar,...
    tank_cap,zeros(sample_q,1),0,E_bSingle,1,Q_Ss,on0,UpPenal,DoPenal,quarters);

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

% Calculate investment costs 
%---------------------------

% Optimal is of course one CHP. Every CHP you add will have a slightly
% higher investment cost (imbalance remaining is less and thus more
% expensive to compensate). This part of the code let you calculate a
% random case (you can select one from the output of the previous part).

Cust = N; % Custopt;

%% RESULTS 
%-----------------------------------
%LOADING RESULTS
%Reading parameters
%-----------------------------------

disp('Loading results day-ahead optimisation...')

if Calc == 1
    wgdx('inputs', i,n,s,P_g,P_c,P_i,P_n,P_st,E_i,E_is,Q_H,Nb,Ns,Ae0,Aq0,Ae1,Aq1,Pi_s,Ecap_lo,Ecap_up,Qcap_up,Cs,bid,bid_bool,bid_single,CHP_bool,QSs,ON0,UpPen,DoPen,Dur,dt);
    gams('CHP'); % Day-ahead optimisation
    da.name = 'day-ahead optimisation';
    [da.obj,da.R_b,da.R_ir,da.FC_bc,da.m_fCHP,da.m_fB,da.Q_CHP,da.Q_B,da.DeltaQ_S,da.Q_S,da.E_CHP,da.E_ir,da.E_b,da.ON] = GAMSREAD(time_i,sample_q,N,S,sample_i);
    try
        TEST(Sample, da, 0, 1)
    catch
        disp('Relax optimisation constraints (on/off penalty; on0)')
        UpPen.val = 0;
        DoPen.val = 0;
        ON0.val = zeros(size(ON0.val));
        wgdx('inputs', i,n,s,P_g,P_c,P_i,P_n,P_st,E_i,E_is,Q_H,Nb,Ns,Ae0,Aq0,Ae1,Aq1,Pi_s,Ecap_lo,Ecap_up,Qcap_up,Cs,bid,bid_bool,bid_single,CHP_bool,QSs,ON0,UpPen,DoPen,Dur,dt);
        gams('CHP'); % Day-ahead optimisation
        da.name = 'day-ahead optimisation';
        [da.obj,da.R_b,da.R_ir,da.FC_bc,da.m_fCHP,da.m_fB,da.Q_CHP,da.Q_B,da.DeltaQ_S,da.Q_S,da.E_CHP,da.E_ir,da.E_b,da.ON] = GAMSREAD(time_i,sample_q,N,S,sample_i);
        try
            TEST(Sample, da, 0, 1)
        catch
            disp('Increase duration')
            Dur.val = 500;
            wgdx('inputs', i,n,s,P_g,P_c,P_i,P_n,P_st,E_i,E_is,Q_H,Nb,Ns,Ae0,Aq0,Ae1,Aq1,Pi_s,Ecap_lo,Ecap_up,Qcap_up,Cs,bid,bid_bool,bid_single,CHP_bool,QSs,ON0,UpPen,DoPen,Dur,dt);
            gams('CHP'); % Day-ahead optimisation
            da.name = 'day-ahead optimisation';
            [da.obj,da.R_b,da.R_ir,da.FC_bc,da.m_fCHP,da.m_fB,da.Q_CHP,da.Q_B,da.DeltaQ_S,da.Q_S,da.E_CHP,da.E_ir,da.E_b,da.ON] = GAMSREAD(time_i,sample_q,N,S,sample_i);
        end
    end
else
    da = resultsOld(nR).da;
end

%-----------------------------------
%POSTPROCESSING RESULTS
%Postproces data
%-----------------------------------

disp('Postprocessing results day-ahead optimisation...')

daT.name = 'actuals (multiple weeks)';
daT.R_b = [daT.R_b; da.R_b(sample_ii)];
daT.R_ir = [daT.R_ir; da.R_ir(sample_ii)];
daT.FC_bc = [daT.FC_bc; da.FC_bc(sample_ii)];
daT.m_fCHP = [daT.m_fCHP; da.m_fCHP(sample_ii,:,:)];
daT.m_fB = [daT.m_fB; da.m_fB(sample_ii,:,:)];
daT.Q_CHP = [daT.Q_CHP; da.Q_CHP(sample_ii,:,:)];
daT.Q_B = [daT.Q_B; da.Q_B(sample_ii,:,:)];
daT.DeltaQ_S = [daT.DeltaQ_S; da.DeltaQ_S(sample_ii,:,:)];
daT.Q_S = [daT.Q_S; da.Q_S(sample_ii,:,:)];
daT.E_CHP = [daT.E_CHP; da.E_CHP(sample_ii,:,:)];
daT.E_ir = [daT.E_ir; da.E_ir(sample_ii,:)];
daT.E_b = [daT.E_b; da.E_b(sample_ii)];
daT.ON = [daT.ON; da.ON(sample_ii,:,:)];

daTT.name = 'actuals (complete simulation)';
daTT.R_b = [daTT.R_b; da.R_b];
daTT.R_ir = [daTT.R_ir; da.R_ir];
daTT.FC_bc = [daTT.FC_bc; da.FC_bc];
daTT.m_fCHP = [daTT.m_fCHP; da.m_fCHP];
daTT.m_fB = [daTT.m_fB; da.m_fB];
daTT.Q_CHP = [daTT.Q_CHP; da.Q_CHP];
daTT.Q_B = [daTT.Q_B; da.Q_B];
daTT.DeltaQ_S = [daTT.DeltaQ_S; da.DeltaQ_S];
daTT.Q_S = [daTT.Q_S; da.Q_S];
daTT.E_CHP = [daTT.E_CHP; da.E_CHP];
daTT.E_ir = [daTT.E_ir; da.E_ir];
daTT.E_b = [daTT.E_b; da.E_b];
daTT.ON = [daTT.ON; da.ON];

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

[Sample] = INITSAMPLE(Par2,Gen2,Year,xi,xj);

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
heatD_s = Sample.heatD_s;

quarters = Gen2.quarters;         
hours = Gen2.hours;               
days_w = Gen2.days_w;             
days_y = Gen2.days_y;             
timestep = Gen2.timestep;         
sample_h = Gen2.sample_h;         
sample_q = Gen2.sample_q;         
time_i = Gen2.time_i;             
time_q = Gen2.time_q;             
tank_cap = Gen2.tank_cap;         
Ecap_upVar = Gen2.Ecap_upVar;     
Ecap_loVar = Gen2.Ecap_loVar;     
Qcap_upVar = Gen2.Qcap_upVar;

if xj ~= 1; Q_Ss = Q_SsDay; end
if xj ~= 1; on0 = on0Day; end

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
[i,n,s,P_g,P_c,P_i,P_n,P_st,E_i,E_is,Q_H,Nb,Ns,Ae0,Aq0,Ae1,Aq1,Pi_s,...
    Ecap_lo,Ecap_up,Qcap_up,Cs,bid,bid_bool,bid_single,CHP_bool,QSs,ON0,UpPen,DoPen,Dur,dt] = GAMSWRITE(sample_q,N,...
    S,price_elecS_s,price_elecC_s,price_actual_s,price_gas_s,actuals_s,...
    heatD_s,Par.Nb0,Par.Ns0,AE0,AQ0,AE1,AQ1,Pi_actual,Ecap_loVar,Ecap_upVar,Qcap_upVar,...
    tank_cap,zeros(sample_q,1),0,Par.E_bSingle,1,Q_Ss,on0,UpPenal,DoPenal,quarters);

BUYst.val=[ones(1,Cust) zeros(1,N-Cust)];
    wgdx('inputs', i,n,s,P_g,P_c,P_i,P_n,P_st,E_i,E_is,Q_H,Nb,Ns,Ae0,Aq0,Ae1,Aq1,Pi_s,Ecap_lo,Ecap_up,Qcap_up,Cs,bid,bid_bool,bid_single,CHP_bool,QSs,ON0,UpPen,DoPen,Dur,dt);

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
[i,n,s,P_g,P_c,P_i,P_n,P_st,E_i,E_is,Q_H,Nb,Ns,Ae0,Aq0,Ae1,Aq1,Pi_s,...
    Ecap_lo,Ecap_up,Qcap_up,Cs,bid,bid_bool,bid_single,CHP_bool,QSs,ON0,UpPen,DoPen,Dur,dt] = GAMSWRITE(sample_q,N,...
    S,price_elecS_s,price_elecC_s,price_actual_s,price_gas_s,actuals_s,...
    heatD_s,Par.Nb0,Par.Ns0,AE0,AQ0,AE1,AQ1,Pi_actual,Ecap_loVar,Ecap_upVar,Qcap_upVar,...
    tank_cap,da.E_b((24-Dnom)*quarters+1:end),1,Par.E_bSingle,1,Q_Ss,on0,UpPenal,DoPenal,quarters);
    

if Calc == 1
    wgdx('inputs', i,n,s,P_g,P_c,P_i,P_n,P_st,E_i,E_is,Q_H,Nb,Ns,Ae0,Aq0,Ae1,Aq1,Pi_s,Ecap_lo,Ecap_up,Qcap_up,Cs,bid,bid_bool,bid_single,CHP_bool,QSs,ON0,UpPen,DoPen,Dur,dt);
    gams('CHP'); % Calculate actuals
    a.name = 'actuals (with CHP; day-ahead bidding)';
    [a.obj,a.R_b,a.R_ir,a.FC_bc,a.m_fCHP,a.m_fB,a.Q_CHP,a.Q_B,a.DeltaQ_S,a.Q_S,a.E_CHP,a.E_ir,a.E_b,a.ON] = GAMSREAD(time_i,sample_q,N,S,sample_i);
    try
        TEST(Sample, a, 0, 1)
    catch
        disp('Relax optimisation constraints (on/off penalty; on0)')
        UpPen.val = 0;
        DoPen.val = 0;
        ON0.val = zeros(size(ON0.val));
        wgdx('inputs', i,n,s,P_g,P_c,P_i,P_n,P_st,E_i,E_is,Q_H,Nb,Ns,Ae0,Aq0,Ae1,Aq1,Pi_s,Ecap_lo,Ecap_up,Qcap_up,Cs,bid,bid_bool,bid_single,CHP_bool,QSs,ON0,UpPen,DoPen,Dur,dt);
        gams('CHP'); % Calculate actuals
        a.name = 'actuals (with CHP; day-ahead bidding)';
        [a.obj,a.R_b,a.R_ir,a.FC_bc,a.m_fCHP,a.m_fB,a.Q_CHP,a.Q_B,a.DeltaQ_S,a.Q_S,a.E_CHP,a.E_ir,a.E_b,a.ON] = GAMSREAD(time_i,sample_q,N,S,sample_i);
        TEST(Sample, a, 0, 1)
    end
else
    a = resultsOld(nR).a;
end

% Actuals (without CHP; day-ahead bidding/optimal)
CHPBool = 0;
[i,n,s,P_g,P_c,P_i,P_n,P_st,E_i,E_is,Q_H,Nb,Ns,Ae0,Aq0,Ae1,Aq1,Pi_s,...
    Ecap_lo,Ecap_up,Qcap_up,Cs,bid,bid_bool,bid_single,CHP_bool,QSs,ON0,UpPen,DoPen,Dur,dt] = GAMSWRITE(sample_q,N,...
    S,price_elecS_s,price_elecC_s,price_actual_s,price_gas_s,actuals_s,...
    heatD_s,Par.Nb0,Par.Ns0,AE0,AQ0,AE1,AQ1,Pi_actual,Ecap_loVar,Ecap_upVar,Qcap_upVar,...
    tank_cap,zeros(sample_q,1),0,Par.E_bSingle,1,Q_Ss,on0,UpPenal,DoPenal,quarters);

BUYst.val=[ones(1,Cust) zeros(1,N-Cust)];
    wgdx('inputs', i,n,s,P_g,P_c,P_i,P_n,P_st,E_i,E_is,Q_H,Nb,Ns,Ae0,Aq0,Ae1,Aq1,Pi_s,Ecap_lo,Ecap_up,Qcap_up,Cs,bid,bid_bool,bid_single,CHP_bool,QSs,ON0,UpPen,DoPen,Dur,dt);

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

SampleAT.sample_i = [SampleAT.sample_i, sample_i(sample_ii)];
SampleAT.imbal_s = [SampleAT.imbal_s; imbal_s(sample_ii,:)];
SampleAT.Pi_st = [SampleAT.Pi_st; Pi_st'];
SampleAT.actuals_s = [SampleAT.actuals_s; actuals_s(sample_ii)];
SampleAT.price_actual_s = [SampleAT.price_actual_s; price_actual_s(sample_ii)];
SampleAT.price_posActual_s = [SampleAT.price_posActual_s; price_posActual_s(sample_ii)];
SampleAT.price_negActual_s = [SampleAT.price_negActual_s; price_negActual_s(sample_ii)];
SampleAT.price_gas_s = [SampleAT.price_gas_s; price_gas_s(sample_ii)];
SampleAT.price_elecC_s = [SampleAT.price_elecC_s; price_elecC_s(sample_ii)];
SampleAT.price_elecS_s = [SampleAT.price_elecS_s; price_elecS_s(sample_ii)];
SampleAT.elecD_s = [SampleAT.elecD_s; elecD_s(sample_ii,:)];
SampleAT.elecWF_s = [SampleAT.elecWF_s; elecWF_s(sample_ii)];
SampleAT.elecWP_s = [SampleAT.elecWP_s; elecWP_s(sample_ii)];
SampleAT.heatD_s = [SampleAT.heatD_s; heatD_s(sample_ii,:)];

SampleATT.sample_i = [SampleATT.sample_i, sample_i];
SampleATT.imbal_s = [SampleATT.imbal_s; imbal_s];
SampleATT.Pi_st = [SampleATT.Pi_st; Pi_st'];
SampleATT.actuals_s = [SampleATT.actuals_s; actuals_s];
SampleATT.price_actual_s = [SampleATT.price_actual_s; price_actual_s];
SampleATT.price_posActual_s = [SampleATT.price_posActual_s; price_posActual_s];
SampleATT.price_negActual_s = [SampleATT.price_negActual_s; price_negActual_s];
SampleATT.price_gas_s = [SampleATT.price_gas_s; price_gas_s];
SampleATT.price_elecC_s = [SampleATT.price_elecC_s; price_elecC_s];
SampleATT.price_elecS_s = [SampleATT.price_elecS_s; price_elecS_s];
SampleATT.elecD_s = [SampleATT.elecD_s; elecD_s];
SampleATT.elecWF_s = [SampleATT.elecWF_s; elecWF_s];
SampleATT.elecWP_s = [SampleATT.elecWP_s; elecWP_s];
SampleATT.heatD_s = [SampleATT.heatD_s; heatD_s];  

aT.name = 'actuals (multiple weeks)';
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

aTT.name = 'actuals (complete simulation)';
aTT.R_b = [aTT.R_b; a.R_b];
aTT.R_ir = [aTT.R_ir; a.R_ir];
aTT.FC_bc = [aTT.FC_bc; a.FC_bc];
aTT.m_fCHP = [aTT.m_fCHP; a.m_fCHP];
aTT.m_fB = [aTT.m_fB; a.m_fB];
aTT.Q_CHP = [aTT.Q_CHP; a.Q_CHP];
aTT.Q_B = [aTT.Q_B; a.Q_B];
aTT.DeltaQ_S = [aTT.DeltaQ_S; a.DeltaQ_S];
aTT.Q_S = [aTT.Q_S; a.Q_S];
aTT.E_CHP = [aTT.E_CHP; a.E_CHP];
aTT.E_ir = [aTT.E_ir; a.E_ir];
aTT.E_b = [aTT.E_b; a.E_b];
aTT.ON = [aTT.ON; a.ON];

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

Q_SsDnom = a.Q_S(Par.Dnom*4,:);
Q_SsDay = a.Q_S(24*4,:);

on0Dnom = a.ON(Par.Dnom*4,:);
on0Day = a.ON(24*4,:);

disp(' ')

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