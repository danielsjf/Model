function [i,n,s,P_g,P_c,P_i,P_n,P_st,E_i,E_is,Q_H,Nb,Ns,Ae0,Aq0,Ae1,Aq1,Pi_s,...
    Ecap_lo,Ecap_up,Qcap_up,Cs,bid,bid_bool,bid_single,CHP_bool,QSs,ON0,UpPen,DoPen,Dur,dt] = GAMSWRITE(sample_q,N,...
    S,price_elecS_s,price_elecC_s,price_imbal_s,price_gas_s,imbal_s,...
    heatD_s,Nb0,Ns0,AE0,AQ0,AE1,AQ1,Pi_st,Ecap_loVar,Ecap_upVar,Qcap_upVar,...
    tank_cap,E_b,E_bBool,E_bSingle,CHPBool,Q_Ss,on0,UpPenal,DoPenal,quarters)
%Input sets
%----------

i.name='i'; %Time
%final.val= 1:15;
i.uels = {{1:sample_q}};             
i.type ='set';

n.name='n'; %Unit
%unit.val= 1:15;
n.uels = {{1:N}};             
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

P_i.name='P_i'; % Imbalance price overproduction
P_i.val= price_imbal_s;
P_i.form = 'full';
P_i.type = 'parameter';
P_i.dim =2;

% P_ineg.name='P_ineg'; % Imbalance price overconsumption
% P_ineg.val= price_negImbal_s;
% P_ineg.form = 'full';
% P_ineg.type = 'parameter';
% P_ineg.dim =1;

P_n.name='P_n'; %price natural gas
P_n.val=price_gas_s;
P_n.form = 'full';
P_n.type = 'parameter';
P_n.dim =1;

P_st.name='P_st'; %price startup
P_st.val=100*ones(N,1);
P_st.form = 'full';
P_st.type = 'parameter';
P_st.dim =1;


E_i.name='E_i'; %electricity demand imbalance
E_i.val=imbal_s;
E_i.form = 'full';
E_i.type = 'parameter';
E_i.dim =2;

E_is.name='E_is'; %electricity demand imbalance
E_is.val=sign(imbal_s);
E_is.form = 'full';
E_is.type = 'parameter';
E_is.dim =2;

Q_H.name='Q_H'; %heat demand house
Q_H.val=heatD_s;
Q_H.form = 'full';
Q_H.type = 'parameter';
Q_H.dim =2;


Nb.name='Nb'; %thermal efficiency boiler
Nb.val=Nb0*ones(N,1);
Nb.form = 'full';
Nb.type = 'parameter';
Nb.dim =1;

Ns.name='Ns'; %thermal efficiency storage
Ns.val=Ns0*ones(N,1);
Ns.form = 'full';
Ns.type = 'parameter';
Ns.dim =1;

Ae0.name='Ae'; %electrical efficiency CHP
Ae0.val=AE0';
Ae0.form = 'full';
Ae0.type = 'parameter';
Ae0.dim =1;

Aq0.name='Aq'; %thermal efficiency CHP
Aq0.val=AQ0';
Aq0.form = 'full';
Aq0.type = 'parameter';
Aq0.dim =1;

Ae1.name='Ae1'; %electrical efficiency CHP
Ae1.val=AE1';
Ae1.form = 'full';
Ae1.type = 'parameter';
Ae1.dim =1;

Aq1.name='Aq1'; %thermal efficiency CHP
Aq1.val=AQ1';
Aq1.form = 'full';
Aq1.type = 'parameter';
Aq1.dim =1;

Pi_s.name='Pi_s'; %probability scenario
Pi_s.val=Pi_st;
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
Cs.val=tank_cap;
Cs.form = 'full';
Cs.type = 'parameter';
Cs.dim =1;

bid.name='bid'; % Bid on the day-ahead market
bid.val=E_b;
bid.form = 'full';
bid.type = 'parameter';
bid.dim=1;

bid_bool.name='bid_bool'; % Take the day-ahead bid into account in the optimalisation (1 is yes, 0 is no)
bid_bool.val=E_bBool;
bid_bool.type = 'parameter';
bid_bool.dim=0;

bid_single.name='bid_single'; % Take the day-ahead bid into account in the optimalisation (1 is yes, 0 is no)
bid_single.val=E_bSingle;
bid_single.type = 'parameter';
bid_single.dim=0;

CHP_bool.name='CHP_bool'; % Take the CHP into account in the optimalisation (1 if CHP on, 0 if CHP off)
CHP_bool.val=CHPBool;
CHP_bool.type = 'parameter';
CHP_bool.dim=0;

QSs.name='Q_Ss'; % tank value in the beginning
QSs.val=Q_Ss;
QSs.form = 'full';
QSs.type = 'parameter';
QSs.dim =1;

MinUp.name='MinUp'; % tank value in the beginning
MinUp.val=4;
MinUp.type = 'parameter';
MinUp.dim =0;

ON0.name='ON0'; % ON value in the beginning
ON0.val=on0;
ON0.form = 'full';
ON0.type = 'parameter';
ON0.dim =1;

UpPen.name='UpPen'; % up value in the beginning
UpPen.val=UpPenal;
UpPen.type = 'parameter';
UpPen.dim =0;

DoPen.name='DoPen'; % dowm value in the beginning
DoPen.val=DoPenal;
DoPen.type = 'parameter';
DoPen.dim =0;

Dur.name='Dur'; % dowm value in the beginning
Dur.val=100;
Dur.type = 'parameter';
Dur.dim =0;

dt.name='dt'; %delta t for the storage time period analyzed 0.25for 15 min
dt.val=1/quarters;
dt.type = 'parameter';
dt.dim=0;

end
