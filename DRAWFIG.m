function [] = DRAWFIG(RES,time_q,scenarios_s,heatD_s,Pi_st,N_sh,S_sh)

if nargin < 5
  error('DRAWFIG:  a, time_q, imbal_s, heatD_s and Pi_st are required inputs')
end

if nargin < 6
  N_sh = 1;
end
N_T = size(heatD_s,2);

if nargin < 7
  S_sh = 1;
end
S_T = size(scenarios_s,2);

GEN = STATSSCEN(RES);
profit = GEN.profit;

E_b = RES.E_b; % Bidding energy demand
FC_bc = RES.FC_bc;
R_b = RES.R_b;
R_ir = RES.R_ir;

imbal_s_S = scenarios_s(:,S_sh); % Wind imbalance (for a certain scenario)
DeltaQ_S_S = RES.DeltaQ_S(:,:,S_sh); % Storage heat supply (for a certain scenario)
Q_S_S = RES.Q_S(:,:,S_sh); % Storage level (for a certain scenario)
Q_B_S = RES.Q_B(:,:,S_sh); % Boiler heat supply (for a certain scenario)
Q_CHP_S = RES.Q_CHP(:,:,S_sh); % CHP heat supply (for a certain scenario)
E_CHP_S = RES.E_CHP(:,:,S_sh); % CHP energy supply (for a certain scenario)
E_i_S = RES.E_i(:,S_sh); % Imbalance reduction (for a certain scenario)

heatD_s_N = heatD_s(:,N_sh); % House heat demand (for a certain unit)

DeltaQ_S_SN = RES.DeltaQ_S(:,N_sh,S_sh); % Storage heat supply (for a certain unit and scenario)
Q_S_SN = RES.Q_S(:,N_sh,S_sh); % Storage level (for a certain unit and scenario)
Q_B_SN = RES.Q_B(:,N_sh,S_sh); % Boiler heat supply (for a certain unit and scenario)
Q_CHP_SN = RES.Q_CHP(:,N_sh,S_sh); % CHP heat supply (for a certain unit and scenario)
E_CHP_NS = RES.E_CHP(:,N_sh,S_sh); % CHP heat supply (for a certain unit and scenario)

m_fB_TS = sum(RES.m_fB(:,:,S_sh),2); % Total fuelflow boilers
m_fCHP_TS = sum(RES.m_fCHP(:,:,S_sh),2); % Total fuelflow CHP's

imbal_s_W = scenarios_s*Pi_st; % Weighted average imbalance
E_i_W = RES.E_i*Pi_st; % Weighted average imbalance reduction

E_CHP_TW = sum(RES.E_CHP,2);
if ~ismatrix(E_CHP_TW)
    E_CHP_TW = permute(E_CHP_TW,[1 3 2]);
end
E_CHP_TW = E_CHP_TW*Pi_st;

figure()
plot(time_q,[imbal_s_S, E_b, E_i_S])
title(['Imbalance for scenario ',num2str(S_sh),' of ',num2str(S_T)]);
legend('Wind imbalance','Bidding','Imbalance reduction');
xlabel('Time [h]');
ylabel('Energy [MWh]');

figure()
subplot(2,1,1)
plot(time_q,Q_S_S)
title(['Storage level for scenario ',num2str(S_sh),' of ',num2str(S_T),' (all units)']);
xlabel('Time [h]');
ylabel('Heat [MWh]');
subplot(2,1,2)
plot(time_q,DeltaQ_S_S)
title(['Storage heat supply for scenario ',num2str(S_sh),' of ',num2str(S_T),' (all units)']);
xlabel('Time [h]');
ylabel('Heat [MWh]');

figure()
plot(time_q,[Q_S_SN,DeltaQ_S_SN]);
title(['Storage unit ',num2str(N_sh),' of ',num2str(N_T),' for scenario ',num2str(S_sh),' of ',num2str(S_T)]);
legend('Storage level','Storage heat supply');
xlabel('Time [h]');
ylabel('Heat [MWh]');

figure()
subplot(2,1,1)
plot(time_q,E_CHP_S)
title(['CHP energy supply for scenario ',num2str(S_sh),' of ',num2str(S_T),' (all units)']);
xlabel('Time [h]');
ylabel('Energy [MWh]');
subplot(2,1,2)
plot(time_q,Q_CHP_S)
title(['CHP heat supply for scenario ',num2str(S_sh),' of ',num2str(S_T),' (all units)']);
xlabel('Time [h]');
ylabel('Heat [MWh]');

figure()
plot(time_q,[heatD_s_N,Q_B_SN,Q_CHP_SN])
title(['Heating balance for unit ',num2str(N_sh),' of ',num2str(N_T),' for scenario ',num2str(S_sh),' of ',num2str(S_T)]);
legend('House heat demand','Boiler heat supply','CHP heat supply');
xlabel('Time [h]');
ylabel('Heat [MWh]');

figure()
plot(time_q,[m_fB_TS,m_fCHP_TS])
title(['Total fuelflow for scenario ',num2str(S_sh),' of ',num2str(S_T)]);
legend('Total fuelflow Boilers','Total fuelflow CHPs');
xlabel('Time [h]');
ylabel('Fuelflow [MWh]');

figure()
plot(time_q,[imbal_s_S,E_CHP_NS,E_b,E_i_S])
title(['Energy balance for unit ',num2str(N_sh),' of ',num2str(N_T),' for scenario ',num2str(S_sh),' of ',num2str(S_T)]);
legend('Imbalance','CHP energy supply unit 1','Bidding','Total imbalance reduction');
xlabel('Time [h]');
ylabel('Energy [MWh]');

figure()
plot(time_q,[imbal_s_W,E_CHP_TW,E_b,E_i_W])
title(['Weighted average (scenarios) total (units) energy balance']);
legend('Weighted average imbalance','CHP total weighted average energy supply','Bidding','Weighted average imbalance reduction');
xlabel('Time [h]');
ylabel('Energy [MWh]');

figure()
plot(time_q,[FC_bc,R_b,R_ir,profit]);
title('Total profit');
legend('Fuel consumption cost','Bidding revenue','Imbalance reduction revenue','profit');
xlabel('Time [h]');
ylabel('Money [€]');

% figure()
% plot(imbal_s,'DisplayName','imbal_s');hold all;plot(actuals_s,'DisplayName','actuals_s');hold off;

end