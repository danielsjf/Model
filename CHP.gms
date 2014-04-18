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
     bid(i)       bid on the day-ahead market
     bid_bool     bid done on the market (0 is false and 1 is true)
     bid_single   single bid for the whole day (0 is multiple and 1 is single)
     CHP_bool     CHP on-off (0 is off and 1 is on)
     BUYst(n)     startvalue CHP's;


Scalar
     dt           timestep;


$gdxin inputs
$load i n s P_g P_c P_i P_n P_st E_i0 E_i1 Q_H Nb Ns Ae Aq Pi_s Ecap_lo Ecap_up Qcap_up Cs bid bid_bool bid_single CHP_bool dt
$gdxin

*$if exist matdata.gms  $include matdata.gms

variables
     obj               objective

     m_fCHP(i,n,s)     primary fuelflow CHP
     m_fB(i,n,s)       primary fuelflow boiler

     Q_CHP(i,n,s)      heat supply CHP
     Q_B(i,n,s)        heat supply boiler
     Q_S(i,n,s)        heat storage
     DeltaQ_S(i,n,s)   heat supply storage

     R_b               Revenue bidding
     R_ir              Revenue imbalance reduction

     C_ef              Cost extra fuel

     FC_b              Fuelcost boiler
     FC_bc             Fuelcost boiler and chp

     E_CHP(i,n,s)      electricity supply CHP
     E_i(i,s)          imbalance reduction
     E_b(i)            electricity demand bidding

     ON(i,n,s)         Turn CHP on
     BUY(n)            Invest in CHP;

     E_i.lo(i,s) = min(0,-E_i0(i,s));
     E_i.up(i,s) = max(-E_i0(i,s),0);

     DeltaQ_S.lo(i,n,s) = -Qcap_up(n);
     DeltaQ_S.up(i,n,s) = Qcap_up(n);


positive variables m_fCHP(i,n,s),m_fB(i,n,s),E_CHP(i,n,s),Q_CHP(i,n,s),Q_B(i,n,s),Q_S(i,n,s),E_b;
     m_fCHP.up(i,n,s) = Ecap_up(n)/Ae(n);
     m_fB.up(i,n,s) = Qcap_up(n)/Nb(n);

     E_CHP.up(i,n,s) = Ecap_up(n);

     Q_CHP.up(i,n,s) = Ecap_up(n)/Ae(n)*Aq(n);
     Q_B.up(i,n,s)= Qcap_up(n);

     E_b.up(i) = sum(n, Ecap_up(n));
*     Q_S.fx('192',n,s)= 0;
     Q_S.up(i,n,s)= Cs(n);

     E_b.fx(i)$(bid_bool) = bid(i);
     E_b.fx(i)$(not CHP_bool) = 0;
     m_fCHP.fx(i,n,s)$(not CHP_bool) = 0


binary variable ON(i,n,s),BUY(n);


*semicont variable E_CHP(i,n),Q_CHP(i,n);
*     E_CHP.lo(i,n) = Ecap_lo(n);
*     E_CHP.up(i,n) = Ecap_up(n);
*
*     Q_CHP.lo(i,n) = Ecap_lo(n)/Ae(n)*Aq(n);
*     Q_CHP.up(i,n) = Ecap_up(n)/Ae(n)*Aq(n);





equation Investment,Revenue_bidding,Revenue_imbalance_red,Cost_extra_fuel,Fuelcost_only_boiler,Fuelcost_boiler_chp,Heat_demand(i,n,s),Grid(i,s),Bidding(i),Storage(i,n,s),CHP_E(i,n,s),CHP_Q(i,n,s),Boiler(i,n,s),Shutdown1(i,n,s),Shutdown2(i,n,s);
*, start_shut1,min_up1,Dispatch(i);

   Investment..                obj =e= sum(i, R_b(i) + R_ir(i) - C_ef(i));

   Revenue_bidding(i)..        R_b(i) =e= E_b(i)*P_g(i);

   Revenue_imbalance_red(i)..  R_ir(i) =e= sum(s, -Pi_s(s)*E_i1(i,s)*E_i(i,s)*P_i(i,s));

   Cost_extra_fuel(i)..        C_ef(i) =e= FC_bc(i);
*- FC_b;

   Fuelcost_only_boiler(i)..   FC_b(i) =e= sum((n,s), Pi_s(s)*(Q_H(i,n)/Nb(n))*P_n(i));

   Fuelcost_boiler_chp(i)..    FC_bc(i) =e= sum((n,s), Pi_s(s)*(m_fCHP(i,n,s)+m_fB(i,n,s))*P_n(i));

   Heat_demand(i,n,s)..        Q_H(i,n) =e= DeltaQ_S(i,n,s) + Q_B(i,n,s) + Q_CHP(i,n,s);

   Grid(i,s)..                 sum(n, E_CHP(i,n,s)) =e= E_i(i,s) + E_b(i);

   Bidding(i)..                E_b(i)$(bid_single) =e= E_b(i--1)$(bid_single);

   Storage(i,n,s)..            Q_S(i,n,s) =e= Q_S(i--1,n,s)*Ns(n) - DeltaQ_S(i,n,s);

   CHP_E(i,n,s)..              E_CHP(i,n,s) =e= m_fCHP(i,n,s)*Ae(n);

   CHP_Q(i,n,s)..              Q_CHP(i,n,s) =e= m_fCHP(i,n,s)*Aq(n);

   Boiler(i,n,s)..             Q_B(i,n,s) =e= m_fB(i,n,s)*Nb(n);

   Shutdown1(i,n,s)..          E_CHP(i,n,s) =g= Ecap_lo(n)*ON(i,n,s);

   Shutdown2(i,n,s)..          E_CHP(i,n,s) =l= Ecap_up(n)*ON(i,n,s);

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
qp1.optcr=0.00001;
qp1.threads=-1;
option MIP=cplex;
*option MIP=gurobi;
solve qp1 using mip maximizing obj;

R_b.l(i)$(not R_b.l(i)) = eps;

R_ir.l(i)$(not R_ir.l(i)) = eps;

FC_bc.l(i)$(not FC_bc.l(i)) = eps;

m_fCHP.l(i,n,s)$(not m_fCHP.l(i,n,s)) = eps;

m_fB.l(i,n,s)$(not m_fB.l(i,n,s)) = eps;

Q_CHP.l(i,n,s)$(not Q_CHP.l(i,n,s)) = eps;

Q_B.l(i,n,s)$(not Q_B.l(i,n,s)) = eps;

Q_S.l(i,n,s)$(not Q_S.l(i,n,s)) = eps;

DeltaQ_S.l(i,n,s)$(not DeltaQ_S.l(i,n,s)) = eps;

E_CHP.l(i,n,s)$(not E_CHP.l(i,n,s)) = eps;

E_i.l(i,s)$(not E_i.l(i,s)) = eps;

E_b.l(i)$(not E_b.l(i)) = eps;



ON.l(i,n,s)$(not ON.l(i,n,s)) = eps;

*BUY.l(n)$(not BUY.l(n)) = eps;

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


execute_unload 'results',  obj, R_b, R_ir, FC_bc, m_fCHP, m_fB, Q_CHP, Q_B, Q_S, DeltaQ_S, E_CHP, E_i, E_b, ON, BUY;
