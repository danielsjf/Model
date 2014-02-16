clear all 
close all
clc
path(path,'C:\GAMS\win64\24.2')
clear gamso;
 load spot2010
 %load ElectPV_multhouse
 load profiles

 
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

%Input constants
%---------------

Cu = 10; % Maximum amount of CHP's
Cust = 9; % Initial amount of CHP's

Cs0 = 200; % [kWh] Storage tank capacity

S = 2; % Number of stochastic scenario's

Ae0 = 0.4; % CHP electrical efficiency
Aq0 = 0.45; % CHP thermal efficiency
Nb0 = 0.98; % Boiler thermal efficiency
Ns0 = 0.998; % Storage heat efficiency

%Input time
%----------

week = 3; % week of the year

quarters = 4; % quarters per hour
hours = 24; % hours per day
days_w = 7; % days per week
days_y = 365; % days per year

timestep = 1/quarters; % [hours] 15 minutes

sample_h = 2*24; % [hours] sample duration
sample_q = sample_h * quarters; % [quarters] sample duration

period_h = week * days_w * hours; % [hours] period of the year
period_q = period_h * quarters; % [quarters] period of the year

%Input time postprocessing
%-------------------------

time_i = (1:sample_q)'; % time index vector sample
time_q = time_i/quarters; % [quarters] time vector sample

sample_i = [1:sample_q] + period_q; % sample period index

total_q=period_q+(sample_q-1); % [quarters] period before and during sample period

%Input price
%-----------

% Gas
gasPrice = 0.040/2; % [€/kWh] gas price (http://epp.eurostat.ec.europa.eu/statistics_explained/index.php/Electricity_and_natural_gas_price_statistics;http://appsso.eurostat.ec.europa.eu/nui/show.do?dataset=nrg_pc_204&lang=en;http://appsso.eurostat.ec.europa.eu/nui/show.do?dataset=nrg_pc_205&lang=en;http://appsso.eurostat.ec.europa.eu/nui/show.do?dataset=nrg_pc_202&lang=en;http://appsso.eurostat.ec.europa.eu/nui/show.do?dataset=nrg_pc_203&lang=en)
price_gas_y = gasPrice*ones(days_y*hours*quarters,1); % [€/kWh] gas price during the year
price_gas_s = price_gas_y(sample_i); % [€/kWh] gas price during the sample duration

% Consumption
price_elecC_d = reshape(repmat([0.15 0.15 0.15 0.15 0.15 0.15 0.15 0.22 0.22 0.22 0.22 0.22 0.22 0.22 0.22 0.22 0.22 0.22 0.22 0.22 0.22 0.22 0.15 0.15],quarters,1),hours*quarters,1); % [€/kWh] electricity price (day/night, one day)
price_elecC_y = repmat(price_elecC_d, days_y, 1); % [€/kWh] consumption electricity price during the year
price_elecC_s = price_elecC_y(sample_i); % [€/kWh] consumption electricity price during the sample duration

% Grid (spotprice)
price_elecS_y = reshape(repmat(spot2010/1000,quarters,1),days_y*hours*quarters,1); % [€/kWh] spot electricity price during the year
price_elecS_s = price_elecS_y(sample_i); % [€/kWh] spot electricity price during the sample duration

% Imbalance
price_imbal_y = repmat(price_elecS_y*2,1,S); % [€/kWh] imbalance electricity price during the year
price_imbal_s = price_imbal_y(sample_i,S); % [€/kWh] imbalance electricity price during the sample duration
 
%Input heat demand
%-----------------

heatD1_y = reshape(repmat(mhouse1h',quarters,1),days_y*hours*quarters,1); % [kWh] heat demand house 1 during the year
heatD1_s = heatD1_y(sample_i); % [kWh] heat demand house 1 during the sample duration

heatD2_y = reshape(repmat(mhouse2h',quarters,1),days_y*hours*quarters,1); % [kWh] heat demand house 2 during the year
heatD2_s = heatD2_y(sample_i); % [kWh] heat demand house 2 during the sample duration

heatD3_y = reshape(repmat(mhouse3h',quarters,1),days_y*hours*quarters,1); % [kWh] heat demand house 3 during the year
heatD3_s = heatD3_y(sample_i); % [kWh] heat demand house 3 during the sample duration

heatD_y = repmat(heatD1_y,1,Cu);
heatD_s = repmat(heatD1_s,1,Cu);

%Input energy demand
%-------------------

% House electricity demand
energy_y = reshape(repmat(mhouse1e+mhouse2e+mhouse3e,quarters,1),days_y*hours*quarters,1);
energy_s = energy_y(sample_i);

% Imbalance
imbal_y = repmat(energy_y,1,S)-mean(energy_y);
imbal_s = imbal_y(sample_i,:);

% Capacity restrictions
Ecap_loVar = 0.2/quarters*ones(Cu,1);
Ecap_upVar = 1.0/quarters*ones(Cu,1);
Qcap_upVar = 100/quarters*ones(Cu,1);

%-----------------------------------
%INITIALIZATION 
%Input variables in structure form 
%-----------------------------------

%Input sets
%----------

i.name='i'; %Time
%final.val= 1:15;
i.uels = {{1:sample_q}};             
i.type ='set';

n.name='n'; %Unit
%unit.val= 1:15;
n.uels = {{1:Cu}};             
n.type ='set';

s.name='s'; %Scenario
%unit.val= 1:15;
s.uels = {{1:S}};             
s.type ='set';

%Input parameters
%----------------

P_g.name='P_g'; %price electricity grid (spotprice)
P_g.val=price_elecS_s; % [€/kWh]
P_g.form = 'full';
P_g.type = 'parameter';
P_g.dim =1;

P_c.name='P_c'; %price electricity consumption
P_c.val= price_elecC_s; % Constant price
P_c.form = 'full';
P_c.type = 'parameter';
P_c.dim =1;

P_i.name='P_i'; %price electricity imbalance
P_i.val= price_imbal_s;
P_i.form = 'full';
P_i.type = 'parameter';
P_i.dim =2;

P_n.name='P_n'; %price natural gas
P_n.val=price_gas_s;
P_n.form = 'full';
P_n.type = 'parameter';
P_n.dim =1;

P_st.name='P_st'; %price startup
P_st.val=100*ones(Cu,1);
P_st.form = 'full';
P_st.type = 'parameter';
P_st.dim =1;


E_i0.name='E_i0'; %electricity demand imbalance
E_i0.val=imbal_s;
E_i0.form = 'full';
E_i0.type = 'parameter';
E_i0.dim =2;

E_i1.name='E_i1'; %electricity demand imbalance
E_i1.val=sign(imbal_s);
E_i1.form = 'full';
E_i1.type = 'parameter';
E_i1.dim =2;

Q_H.name='Q_H'; %heat demand house
Q_H.val=heatD_s;
Q_H.form = 'full';
Q_H.type = 'parameter';
Q_H.dim =2;


Nb.name='Nb'; %thermal efficiency boiler
Nb.val=Nb0*ones(Cu,1);
Nb.form = 'full';
Nb.type = 'parameter';
Nb.dim =1;

Ns.name='Ns'; %thermal efficiency storage
Ns.val=Ns0*ones(Cu,1);
Ns.form = 'full';
Ns.type = 'parameter';
Ns.dim =1;

Ae.name='Ae'; %electrical efficiency CHP
Ae.val=Ae0*ones(Cu,1);
Ae.form = 'full';
Ae.type = 'parameter';
Ae.dim =1;

Aq.name='Aq'; %thermal efficiency CHP
Aq.val=Aq0*ones(Cu,1);
Aq.form = 'full';
Aq.type = 'parameter';
Aq.dim =1;

Pi_s.name='Pi_s'; %probability scenario
Pi_s.val=1/S*ones(S,1);
Pi_s.form = 'full';
Pi_s.type = 'parameter';
Pi_s.dim =1;

Ecap_lo.name='Ecap_lo'; %minimal electric energy supply CHP
Ecap_lo.val=Ecap_loVar;
Ecap_lo.form = 'full';
Ecap_lo.type = 'parameter';
Ecap_lo.dim =1;

Ecap_up.name='Ecap_up'; %maximal electric energy supply CHP
Ecap_up.val=Ecap_upVar;
Ecap_up.form = 'full';
Ecap_up.type = 'parameter';
Ecap_up.dim =1;

Qcap_up.name='Qcap_up'; %maximal thermal energy supply boiler
Qcap_up.val=Qcap_upVar;
Qcap_up.form = 'full';
Qcap_up.type = 'parameter';
Qcap_up.dim =1;

Cs.name='Cs'; %storage tank capacity
Cs.val=Cs0*ones(Cu,1);
Cs.form = 'full';
Cs.type = 'parameter';
Cs.dim =1;

dt.name='dt'; %delta t for the storage time period analyzed 0.25for 15 min
dt.val=1/quarters;
dt.type = 'parameter';
dt.dim=0;

BUYst.name='BUYst'; %storage tank capacity
BUYst.val=[ones(1,Cust) zeros(1,Cu-Cust)];
BUYst.form = 'full';
BUYst.type = 'parameter';
BUYst.dim =1;

index=0;
indexp=0;
index_m=0;

%% CALL GAMS 

%-----------------------------------
%OPTIMIZATION
%Calls GAMS routine 
%-----------------------------------

% Calculate investment costs (optimal)
%-------------------------------------

o2 = 0;
down = 1;

while true
    BUYst.val=[ones(1,Cust) zeros(1,Cu-Cust)];
    wgdx('inputs', i,n,s,P_g,P_c,P_i,P_n,P_st,E_i0,E_i1,Q_H,Nb,Ns,Ae,Aq,Pi_s,Ecap_lo,Ecap_up,Qcap_up,Cs,dt,BUYst);
    path(path,'C:\GAMS\win64\24.2')
    gams('CHP');%LOCAL PRICE IS TAKEN INTO ACCOUNT 
    
    rs.name = 'obj';
    r = rgdx ('results', rs);
    obj=r.val(:,1);
    o1 = obj/Cust;
    [obj Cust o1]
    
    if o1 >= o2
        o2 = o1;
        if down == 1
            Cust = Cust - 1;
            if Cust == 0; break; end;
        else
            Cust = Cust + 1;
            if Cust > Cu; break; end;
        end
    else
        if down == 1
            Cust = Cust + 2;
            if Cust > Cu; break; end;
            down = 0;
        else
            break;
        end
    end
end

Custopt = Cust;

% Calculate investment costs 
%---------------------------

% Optimal is of course one CHP. Every CHP you add will have a slightly
% higher investment cost (imbalance remaining is less and thus more
% expensive to compensate). This part of the code let you calculate a
% random case (you can select one from the output of the previous part).

Cust = Cu; % Custopt;

BUYst.val=[ones(1,Cust) zeros(1,Cu-Cust)];
    wgdx('inputs', i,n,s,P_g,P_c,P_i,P_n,P_st,E_i0,E_i1,Q_H,Nb,Ns,Ae,Aq,Pi_s,Ecap_lo,Ecap_up,Qcap_up,Cs,dt,BUYst);

    gams('CHP');%LOCAL PRICE IS TAKEN INTO ACCOUNT 
    
    rs.name = 'obj';
    r = rgdx ('results', rs);
    obj=r.val(:,1);
    o1 = obj/Cust;
    [obj Cust o1]

%% RESULTS 

%-----------------------------------
%LOADING RESULTS
%Reading parameters
%-----------------------------------

rs.name = 'obj';
r = rgdx ('results', rs);
obj=r.val(:,1); 

rs.name = 'm_fCHP';
r = rgdx ('results', rs);
m_fCHP = zeros(sample_q,Cu);
for j=1:Cu
    m_fCHP(time_i,j)=r.val((time_i-1)*Cu+j,3);
end

rs.name = 'm_fB';
r = rgdx ('results', rs);
m_fB = zeros(sample_q,Cu);
for j=1:Cu
    m_fB(time_i,j)=r.val((time_i-1)*Cu+j,3);
end


rs.name = 'Q_CHP'; 
r = rgdx ('results', rs);
Q_CHP = zeros(sample_q,Cu);
for j=1:Cu
    Q_CHP(time_i,j)=r.val((time_i-1)*Cu+j,3);
end


rs.name = 'Q_B';
r = rgdx ('results', rs);
Q_B = zeros(sample_q,Cu);
for j=1:Cu
    Q_B(time_i,j)=r.val((time_i-1)*Cu+j,3);
end
 
rs.name = 'DeltaQ_S'; %Storage
r = rgdx ('results', rs); 
DeltaQ_S = zeros(sample_q,Cu);
for j=1:Cu
    DeltaQ_S(time_i,j)=r.val((time_i-1)*Cu+j,3);
end

rs.name = 'Q_S';
r = rgdx ('results', rs);
Q_S = zeros(sample_q,Cu);
for j=1:Cu
    Q_S(time_i,j)=r.val((time_i-1)*Cu+j,3);
end


rs.name = 'E_CHP';
r = rgdx ('results', rs);
E_CHP = zeros(sample_q,Cu);
for j=1:Cu
    E_CHP(time_i,j)=r.val((time_i-1)*Cu+j,3);
end

rs.name = 'E_i';
r = rgdx ('results', rs);
E_i = zeros(sample_q,S);
for j=1:S
    E_i(time_i,j)=r.val((time_i-1)*S+j,3);
end

rs.name = 'E_b';
r = rgdx ('results', rs);
E_b=r.val*ones(size(sample_i'));

rs.name = 'ON';
r = rgdx ('results', rs);
ON = zeros(sample_q,Cu);
for j=1:Cu
    ON(time_i,j)=r.val((time_i-1)*Cu+j,3);
end

rs.name = 'BUY';
r = rgdx ('results', rs);
BUY=r.val(:,2);

%-----------------------------------
%POSTPROCESSING RESULTS
%Postproces data
%-----------------------------------

C_f = sum(m_fCHP,2).*P_n.val;

R_b = E_b.*P_g.val;
R_i = abs(E_i*Pi_s.val).*P_i.val;

profit = R_b + R_i - C_f;

profit_t = sum(profit,1) % [€] total profit

%-----------------------------------
%SHOWING RESULTS
%Making figures
%-----------------------------------

figure(1)
plot(time_q,[E_i0.val(:,1), E_b, E_i*Pi_s.val])
title('Imbalance');
legend('Wind imbalance','Bidding','Imbalance reduction');

figure(2)
subplot(2,1,1)
plot(time_q,Q_S)
title('Storage');
subplot(2,1,2)
plot(time_q,DeltaQ_S)
title('Storage heat supply');

figure(3)
plot(time_q,[Q_S(:,1),DeltaQ_S(:,1)]);
legend('Storage','Storage heat supply');
title('Storage unit 1');
xlabel('Time [h]');
ylabel('Heat [kWh]');

figure(4)
subplot(2,1,1)
plot(time_q,E_CHP)
title('CHP energy supply');
xlabel('Time [h]');
ylabel('Energy [kWh]');
subplot(2,1,2)
plot(time_q,Q_CHP)
title('CHP heat supply');
xlabel('Time [h]');
ylabel('Heat [kWh]');

figure(5)
plot(time_q,[Q_H.val(:,1),Q_B(:,1),Q_CHP(:,1)])
legend('House heat demand','Boiler heat supply','CHP heat supply');
title('Heating unit 1');
xlabel('Time [h]');
ylabel('Heat [kWh]');

figure(6)
plot(time_q,[sum(m_fB,2),sum(m_fCHP,2)])
legend('Total fuelflow Boilers','Total fuelflow CHPs');
title('Total fuelflow');
xlabel('Time [h]');
ylabel('Fuelflow [kWh]');

figure(7)
plot(time_q,[E_i0.val(:,1),E_CHP(:,1),E_b(:,1),E_i(:,1)])
legend('Imbalance','CHP energy supply unit 1','Bidding','Total imbalance reduction');
title('Energy unit 1');
xlabel('Time [h]');
ylabel('Energy [kWh]');

figure(8)
plot(time_q,[E_i0.val(:,1),sum(E_CHP,2),E_b(:,1),E_i*Pi_s.val])
legend('Imbalance','CHP energy supply','Bidding','Imbalance reduction');
title('Energy');
xlabel('Time [h]');
ylabel('Energy [kWh]');

figure(9)
plot(time_q,[C_f,R_b,R_i,profit]);
legend('Fuel consumption cost','Bidding revenue','Imbalance reduction revenue','profit');
title('Profit');
xlabel('Time [h]');
ylabel('Money [€]');

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