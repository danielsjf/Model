function [GEN, CHP] = STATSSCEN(RES,PAR,S,dib,Period)
if nargin < 3
  S = 1;
end
if nargin < 4
  dib = 0;
end
if nargin < 2
  error('Not enough arguments');
end
if nargin < 5
  Period = 1:size(RES.E_CHP,1);
end
if S > size(RES.E_CHP,3)
    S = 1;
end

% Compare profit, storage use, boiler use, CHP use (on/off time, on/off
% cycles) for the following scenarios:

if dib == 1, disp(['Data for the ',RES.name]); end
if dib == 1, disp('  '); end

% General
%--------

N = size(RES.E_CHP,2); % number of units
time_l = size(RES.E_CHP(Period,:),1); % length of the timeperiod

% Profit
%-------
GEN.costWchp = RES.R_b(Period) + RES.R_ir(Period) - RES.FC_bc(Period); % [€] cost at every moment in time
GEN.costWchpT = sum(GEN.costWchp,1); % [€] total cost

if dib == 1, disp(['Cost with CHP: ',num2str(GEN.costWchpT)]); end
if dib == 1, disp('  '); end

% CHP
%----
% Average uptime (all units), maximum uptime (all units), minimum uptime
% (all units)

% Binairy matrix (-1 is off, 1 is on)
CHPbin = RES.E_CHP(Period,:,S);
CHPbin(CHPbin>0) = 1;
CHPbin(CHPbin==0) = -1;

% Usage
CHP.used = zeros(1,N); % [%] On-time of the CHP

for k = 1:N
    CHP.used(k) = numel(find(CHPbin(:,k)>0))/time_l*100;
end

CHP.usedM = [find(CHP.used==max(CHP.used)),round(CHP.used(CHP.used==max(CHP.used))*100)/100]; % [unit, %] Most used CHP: unit number and percentage
CHP.usedM = [CHP.usedM(1),CHP.usedM(end)];
CHP.usedL = [find(CHP.used==min(CHP.used)),round(CHP.used(CHP.used==min(CHP.used))*100)/100]; % [unit, %] Least used CHP: unit number and percentage
CHP.usedL = [CHP.usedL(1),CHP.usedL(end)];
CHP.usedA = round(sum(CHP.used)/N*100)/100; % [%] Average usage of a CHP
if dib == 1, disp(['The most used CHP is unit ',num2str(CHP.usedM(1)),' which is on ',num2str(CHP.usedM(2)),'% of the time.']); end
if dib == 1, disp(['The least used CHP is unit ',num2str(CHP.usedL(1)),' which is on ',num2str(CHP.usedL(2)),'% of the time.']); end
if dib == 1, disp(['The average CHP is used ',num2str(CHP.usedA(1)),'% of the time.']); end
if dib == 1, disp('  '); end

% On/off cycles
CHPcum = cumsum(CHPbin,1);
CHP.switches = zeros(1,N); % [#] Number of on and off switches of the CHP
for k = 1:N
    CHP.switches(1,k) = numel(findpeaks(CHPcum(:,k),'THRESHOLD',0))*2 + (CHPbin(2,k)==-1) + (CHPbin(1,k)~=CHPbin(2,k)) + (CHPbin(end,k)==1);
end

CHP.switchesM = [find(CHP.switches==max(CHP.switches)),CHP.switches(CHP.switches==max(CHP.switches))]; % [unit, %] CHP with the most on/off switches: unit number and percentage
CHP.switchesM = [CHP.switchesM(1),CHP.switchesM(end)];
CHP.switchesL = [find(CHP.switches==min(CHP.switches)),CHP.switches(CHP.switches==min(CHP.switches))]; % [unit, %] CHP with the least on/off switches: unit number and percentage
CHP.switchesL = [CHP.switchesL(1),CHP.switchesL(end)];
CHP.switchesA = round(sum(CHP.switches)/N*100)/100; % [%] Average on/off switches of a CHP
if dib == 1, disp(['The CHP which is switched on or off the most is unit ',num2str(CHP.switchesM(1)),' which is switched ',num2str(CHP.switchesM(2)),' times on or off.']); end
if dib == 1, disp(['The CHP which is switched on or off the least is unit ',num2str(CHP.switchesL(1)),' which is switched ',num2str(CHP.switchesL(2)),' times on or off.']); end
if dib == 1, disp(['The average CHP is switched on or off ',num2str(CHP.switchesA),' times.']); end
if dib == 1, disp('  '); end

% Uptime/downtime
CHP.uptimes = zeros(time_l,N);
CHP.uptimesL = zeros(1,N);
CHP.uptimesM = zeros(1,N);
CHP.uptimesA = zeros(1,N);
CHP.downtimes = zeros(time_l,N);
CHP.downtimesL = zeros(1,N);
CHP.downtimesM = zeros(1,N);
CHP.downtimesA = zeros(1,N);
for k = 1:N
    temp = [(CHPbin(1,k)==1)*2-(CHPbin(1,k)==-1)*2;diff(CHPbin(:,k));(CHPbin(end,k)==-1)*2-(CHPbin(end,k)==1)*2]; % Array with local extrema as nonzero elements (-2 for maximum and 2 for minimum)
    temp1 = find(temp); % Elements where a local extrema is reached
    temp2 = diff(temp1); % Distances between local extrema
    if numel(temp2) == 1
        updown = CHPbin(1,k);
        CHP.uptimes(1,k) = 0;
        CHP.uptimesM(k) = 0;
        CHP.uptimesA(k) = 0;
        CHP.downtimes(1,k) = 0;
        CHP.downtimesM(k) = 0;
        CHP.downtimesA(k) = 0;
        if updown == 1
            CHP.uptimes(1,k) = time_l/4;
            CHP.uptimesM(k) = CHP.uptimes(1,k);
            CHP.uptimesA(k) = CHP.uptimes(1,k);
        else
            CHP.downtimes(1,k) = time_l/4;
            CHP.downtimesM(k) = CHP.downtimes(1,k);
            CHP.downtimesA(k) = CHP.downtimes(1,k);
        end
    else
        i = floor(numel(temp2)/2);
        CHP.uptimes(1:floor(numel(temp2)/2)+(CHPbin(1,k)==1&&numel(temp2)/2~=i),k) = temp2(1+(CHPbin(1,k)==-1):2:numel(temp2))/4; % [hours] Distances between local extrema (uptime)
        CHP.uptimesL(k) = min(CHP.uptimes(CHP.uptimes(:,k)>0,k)); % [hours] Smallest distance between local extrema (uptime)
        CHP.uptimesM(k) = max(CHP.uptimes(:,k)); % [hours] Largest distance between local extrema (uptime)
        CHP.uptimesA(k) = mean(CHP.uptimes(CHP.uptimes(:,k)>0)); % [hours] Average distance between local extrema (uptime)
        CHP.downtimes(1:floor(numel(temp2)/2)+(CHPbin(1,k)==-1&&numel(temp2)/2~=i),k) = temp2(1+(CHPbin(1,k)==1):2:numel(temp2))/4; % [hours] Distances between local extrema (downtime)
        CHP.downtimesL(k) = min(CHP.downtimes(CHP.downtimes(:,k)>0,k)); % [hours] Smallest distance between local extrema (downtime)
        CHP.downtimesM(k) = max(CHP.downtimes(:,k)); % [hours] Largest distance between local extrema (downtime)
        CHP.downtimesA(k) = mean(CHP.downtimes(CHP.downtimes(:,k)>0)); % [hours] Average distance between local extrema (downtime)
    end
end

CHP.uptimeL = [find(CHP.uptimesL==min(CHP.uptimesL)),CHP.uptimesL(CHP.uptimesL==min(CHP.uptimesL))]; % [unit, hours] CHP with the shortest uptime: unit number and duration
CHP.uptimeL = [CHP.uptimeL(1),CHP.uptimeL(end)];
CHP.uptimeM = [find(CHP.uptimesM==max(CHP.uptimesM)),CHP.uptimesM(CHP.uptimesM==max(CHP.uptimesM))]; % [unit, hours] CHP with the longest uptime: unit number and duration
CHP.uptimeM = [CHP.uptimeM(1),CHP.uptimeM(end)];
CHP.uptimeA = sum(CHP.uptimesA)/N; % [hours] Average CHP uptime
CHP.downtimeL = [find(CHP.downtimesL==min(CHP.downtimesL)),CHP.downtimesL(CHP.downtimesL==min(CHP.downtimesL))]; % [unit, hours] CHP with the longest downtime: unit number and duration
CHP.downtimeL = [CHP.downtimeL(1),CHP.downtimeL(end)];
CHP.downtimeM = [find(CHP.downtimesM==max(CHP.downtimesM)),CHP.downtimesM(CHP.downtimesM==max(CHP.downtimesM))]; % [unit, hours] CHP with the longest downtime: unit number and duration
CHP.downtimeM = [CHP.downtimeM(1),CHP.downtimeM(end)];
CHP.downtimeA = sum(CHP.downtimesA)/N; % [hours] Average CHP downtime
if dib == 1, disp(['The CHP with the shortest uptime is unit ',num2str(CHP.uptimeL(1)),' which is up for ',num2str(CHP.uptimeL(2)),' hours.']); end
if dib == 1, disp(['The CHP with the longest uptime is unit ',num2str(CHP.uptimeM(1)),' which is up for ',num2str(CHP.uptimeM(2)),' hours.']); end
if dib == 1, disp(['The CHP with the shortest downtime is unit ',num2str(CHP.downtimeL(1)),' which is down for ',num2str(CHP.downtimeL(2)),' hours.']); end
if dib == 1, disp(['The CHP with the longest downtime is unit ',num2str(CHP.downtimeM(1)),' which is down for ',num2str(CHP.downtimeM(2)),' hours.']); end
if dib == 1, disp(['The average CHP uptime is ',num2str(CHP.uptimeA),' hours and the average CHP downtime is ',num2str(CHP.downtimeA),' hours.']); end
if dib == 1, disp('  '); end

% Efficiencies

FuelEffCHP = sum(sum(RES.Q_CHP(Period,:) + RES.E_CHP(Period,:),2))/sum(sum(RES.m_fCHP,2));
FuelEffB = sum(sum(RES.Q_B(Period,:),2))/sum(sum(RES.m_fB(Period,:),2));

StorageLosses = sum(RES.Q_S(Period(1),:)) - sum(sum(RES.DeltaQ_S(Period,:),2)) - sum(RES.Q_S(Period(end),:)); % Storage losses
TotalProduction = sum(sum(RES.Q_CHP(Period,:) + RES.E_CHP(Period,:) + RES.Q_B(Period,:),2)); % Total energy production
FuelEffWOL = (TotalProduction - StorageLosses)/sum(sum(RES.m_fCHP(Period,:) + RES.m_fB(Period,:),2)); % Total fuel efficiency without storage losses
FuelEffWL = TotalProduction/sum(sum(RES.m_fCHP(Period,:) + RES.m_fB(Period,:),2)); % Total fuel efficiency with storage losses
if dib == 1, disp(['The fuel utilisation ratio of the CHP is ',num2str(FuelEffCHP*100,4),'%']); end
if dib == 1, disp(['The fuel utilisation ratio of the Boiler is ',num2str(FuelEffB*100,4),'%']); end
if dib == 1, disp(['The total fuel utilisation ratio without storage losses is ',num2str(FuelEffWOL*100,4),'%']); end
if dib == 1, disp(['The storage losses are ',num2str(StorageLosses*1000,4),' [kWh] on a total production of ',num2str(TotalProduction,4),' [MWh] (heat and electricity)']); end
if dib == 1, disp(['The total fuel utilisation ratio with storage losses is ',num2str(FuelEffWL*100,4),'%']); end

if dib == 1, disp('  '); end

HeatFact = 1 - (15+273)/(70+273);
ExEffCHP = sum((PAR.AE0+PAR.AQ0*HeatFact).*sum(RES.m_fCHP(Period,:),1),2)/sum(sum(RES.m_fCHP(Period,:),2));
ExEffB = sum((PAR.Nb0*HeatFact).*sum(RES.m_fB(Period,:),1),2)/sum(sum(RES.m_fB(Period,:),2));

ExEffWOL = (TotalProduction - StorageLosses)/sum(sum(RES.m_fCHP(Period,:) + RES.m_fB(Period,:),2)); % Total fuel efficiency without storage losses
ExEffWL = (ExEffCHP*sum(sum(RES.m_fCHP(Period,:),2)) + ExEffB*sum(sum(RES.m_fB(Period,:)),2)) / sum(sum(RES.m_fCHP(Period,:) + RES.m_fB(Period,:),2)); % Total fuel efficiency with storage losses
if dib == 1, disp(['The exergetic efficiency of the CHP is ',num2str(ExEffCHP*100,4),'%']); end
if dib == 1, disp(['The exergetic efficiency of the Boiler is ',num2str(ExEffB*100,4),'%']); end
if dib == 1, disp(['The total exergetic efficiency is ',num2str(ExEffWL*100,4),'%']); end

% Boiler
%-------

if dib == 1, disp('  '); end


end