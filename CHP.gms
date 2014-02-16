$ontext
The program inted to minimize the cost of the energy  used by a household.
 Which have installed a CHP unit and back up boiler a
$offtext
*$set matout "'matchp.gdx', Y, X, CHP_electrical1, CHP_electrical2, boiler1, boiler2, qst1, qst2, CHP_primary1, CHP_primary2, Echp_house, Echp_grid, charge1,  charge2 ";
*$set matout "'matchp.gdx', ON1, ON2, ON3, CHP_electrical1, CHP_electrical2, CHP_electrical3, boiler1, boiler2, boiler3, qst1, qst2, qst3, CHP_primary1, CHP_primary2, CHP_primary3, charge1, charge2, charge3, obj ";

$onempty
*$include matglobs.gms


set
     i   time
     n   house or unit
     s   scenario;

parameters
     P_g(i)       price electricity grid (spotprice)
     P_c(i)       price electricity consumption
     P_i(i,s)     price electricity imbalance
     P_n(i)       price natural gas
     P_st(n)      price startup

     E_i0(i,s)    electricity demand imbalance
     E_i1(i,s)    electricity demand imbalance sign

     Q_H(i,n)     heat demand house

     Nb(n)        thermal efficiency boiler
     Ns(n)        thermal efficiency storage
     Ae(n)        electrical efficiency CHP
     Aq(n)        thermal efficiency CHP

     Pi_s(s)      probability scenario

     Ecap_lo(n)   minimal electric energy supply CHP
     Ecap_up(n)   maximal electric energy supply CHP

     Qcap_up(n)   maximal thermal energy supply boiler

     Cs(n)        storage tank capacity
     BUYst(n)     startvalue CHP's;


Scalar
     dt           timestep;


$gdxin inputs
$load i n s P_g P_c P_i P_n P_st E_i0 E_i1 Q_H Nb Ns Ae Aq Pi_s Ecap_lo Ecap_up Qcap_up Cs dt BUYst
$gdxin

*$if exist matdata.gms  $include matdata.gms

variables
     obj               objective

     m_fCHP(i,n)       primary fuelflow CHP
     m_fB(i,n)         primary fuelflow boiler

     Q_CHP(i,n)        heat supply CHP
     Q_B(i,n)          heat supply boiler
     Q_S(i,n)          heat storage
     DeltaQ_S(i,n)     heat supply storage

     E_CHP(i,n)        electricity supply CHP
     E_i(i,s)          imbalance reduction
     E_b               electricity demand bidding

     ON(i,n)           Turn CHP on
     BUY(n)            Invest in CHP;

     E_i.lo(i,s) = min(0,-E_i0(i,s));
     E_i.up(i,s) = max(-E_i0(i,s),0);

     DeltaQ_S.lo(i,n) = -Qcap_up(n);
     DeltaQ_S.up(i,n) = Qcap_up(n);


positive variables m_fCHP(i,n),m_fB(i,n),E_CHP(i,n),Q_CHP(i,n),Q_B(i,n),Q_S(i,n),E_b;
     m_fCHP.up(i,n) = Ecap_up(n)/Ae(n);
     m_fB.up(i,n) = Qcap_up(n)/Nb(n);

     E_CHP.up(i,n) = Ecap_up(n);

     Q_CHP.up(i,n) = Ecap_up(n)/Ae(n)*Aq(n);
     Q_B.up(i,n)= Qcap_up(n);

     E_b.up = sum(n, Ecap_up(n));
     Q_S.fx('192',n)= 0;
     Q_S.up(i,n)= Cs(n);


binary variable ON(i,n),BUY(n);
     BUY.l(n) = BUYst(n);
     ON.fx(i,n)$(BUYst(n) < 0.5) = 0;


*semicont variable E_CHP(i,n),Q_CHP(i,n);
*     E_CHP.lo(i,n) = Ecap_lo(n);
*     E_CHP.up(i,n) = Ecap_up(n);
*
*     Q_CHP.lo(i,n) = Ecap_lo(n)/Ae(n)*Aq(n);
*     Q_CHP.up(i,n) = Ecap_up(n)/Ae(n)*Aq(n);





equation Investment,Heat_demand(i,n),Grid(i),Storage(i,n),Boiler(i,n),CHP_E(i,n),CHP_Q(i,n),Shutdown1(i,n),Shutdown2(i,n);
*, start_shut1,min_up1,Dispatch(i);

Investment..                   obj =e= sum((i,n), (Q_B(i,n)+Q_CHP(i,n))/Nb(n)*P_n(i)) + sum(i, E_b*P_g(i)) + sum((i,s), -Pi_s(s)*E_i1(i,s)*E_i(i,s)*P_i(i,s)) - sum((i,n), (m_fCHP(i,n)+m_fB(i,n))*P_n(i));

   Heat_demand(i,n)..          Q_H(i,n) =e= DeltaQ_S(i,n) + Q_B(i,n) + Q_CHP(i,n);

   Grid(i)..                   sum(n, E_CHP(i,n)) =e= sum(s, Pi_s(s)*E_i(i,s)) + E_b;

   Storage(i,n)..              Q_S(i,n) =e= Q_S(i--1,n)*Ns(n) - DeltaQ_S(i,n);

   CHP_E(i,n)..                E_CHP(i,n) =e= m_fCHP(i,n)*Ae(n);

   CHP_Q(i,n)..                Q_CHP(i,n) =e= m_fCHP(i,n)*Aq(n);

   Boiler(i,n)..               Q_B(i,n) =e= m_fB(i,n)*Nb(n);

   Shutdown1(i,n)..            E_CHP(i,n) =g= Ecap_lo(n)*ON(i,n);

   Shutdown2(i,n)..            E_CHP(i,n) =l= Ecap_up(n)*ON(i,n);

*MINIMUM UP TIME

*start_shut1(i)..             ON1(i)- ON1(i-1)=e= Uup1(i)- Udown1(i);
*min_up1(i)..                 ON1(i+1)=g=Uup1(i);


option limrow=4, limcol=4;
OPTION RESLIM = 200000000;

*display DeltaQ_S, P_c, E_i, Cs
*model qp1 /cost, demand,storage,limit,chp_limit,Electrdem;
model qp1 /all /;
*qp1.Workspace = 30;
*qp1.optfile=1;
option MIP=cplex;
solve qp1 using mip maximizing obj;


m_fCHP.l(i,n)$(not m_fCHP.l(i,n)) = eps;

m_fB.l(i,n)$(not m_fB.l(i,n)) = eps;

Q_CHP.l(i,n)$(not Q_CHP.l(i,n)) = eps;

Q_B.l(i,n)$(not Q_B.l(i,n)) = eps;

Q_S.l(i,n)$(not Q_S.l(i,n)) = eps;

DeltaQ_S.l(i,n)$(not DeltaQ_S.l(i,n)) = eps;

E_CHP.l(i,n)$(not E_CHP.l(i,n)) = eps;

E_i.l(i,s)$(not E_i.l(i,s)) = eps;

E_b.l$(not E_b.l) = eps;

ON.l(i,n)$(not ON.l(i,n)) = eps;

BUY.l(n)$(not BUY.l(n)) = eps;

$ontext
CHP_electrical1.l(i)$(not CHP_electrical1.l(i)) = eps;

boiler1.l(i)$(not boiler1.l(i)) = eps;


qst1.l(i)$(not qst1.l(i)) = eps;


ON1.l(i)$(not ON1.l(i)) = eps;


CHP_primary1.l(i)$(not CHP_primary1.l(i)) = eps;

Echp_grid.l(i)$(not Echp_grid.l(i)) = eps;
Echp_house.l(i)$(not Echp_house.l(i)) = eps;

charge1.l(i)$(not charge1.l(i)) = eps;
$offtext


*execute_unload %matout%;


execute_unload 'results',  obj, m_fCHP, m_fB, Q_CHP, Q_B, Q_S, DeltaQ_S, E_CHP, E_i, E_b, ON, BUY;
