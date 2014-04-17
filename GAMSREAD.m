function [obj,R_b,R_ir,FC_bc,m_fCHP,m_fB,Q_CHP,Q_B,DeltaQ_S,Q_S,E_CHP,E_i,E_b,ON] = GAMSREAD(time_i,sample_q,N,S,sample_i)
% Objective
rs.name = 'obj';
r = rgdx ('results', rs);
obj=r.val(:,1); 

% Bidding revenue
rs.name = 'R_b';
r = rgdx ('results', rs);
R_b(time_i,1)=r.val(time_i,2);

% Imbalance reduction revenue
rs.name = 'R_ir';
r = rgdx ('results', rs);
R_ir(time_i,1)=r.val(time_i,2);

% Boiler and CHP fuelcost
rs.name = 'FC_bc';
r = rgdx ('results', rs);
FC_bc(time_i,1)=r.val(time_i,2);

% CHP fuelflow
rs.name = 'm_fCHP';
r = rgdx ('results', rs);
m_fCHP = zeros(sample_q,N,S);
for j=1:N
    for k=1:S
        m_fCHP(time_i,j,k)=r.val((time_i-1)*N*S+(j-1)*S+k,4);
    end
end

% Boiler fuelflow
rs.name = 'm_fB';
r = rgdx ('results', rs);
m_fB = zeros(sample_q,N,S);
for j=1:N
    for k=1:S
        m_fB(time_i,j,k)=r.val((time_i-1)*N*S+(j-1)*S+k,4);
    end
end

% CHP heat supply
rs.name = 'Q_CHP'; 
r = rgdx ('results', rs);
Q_CHP = zeros(sample_q,N,S);
for j=1:N
    for k=1:S
        Q_CHP(time_i,j,k)=r.val((time_i-1)*N*S+(j-1)*S+k,4);
    end
end

% Boiler heat supply
rs.name = 'Q_B';
r = rgdx ('results', rs);
Q_B = zeros(sample_q,N,S);
for j=1:N
    for k=1:S
        Q_B(time_i,j,k)=r.val((time_i-1)*N*S+(j-1)*S+k,4);
    end
end
 
% Storage heat supply
rs.name = 'DeltaQ_S';
r = rgdx ('results', rs); 
DeltaQ_S = zeros(sample_q,N,S);
for j=1:N
    for k=1:S
        DeltaQ_S(time_i,j,k)=r.val((time_i-1)*N*S+(j-1)*S+k,4);
    end
end

% Storage heat volume
rs.name = 'Q_S';
r = rgdx ('results', rs);
Q_S = zeros(sample_q,N,S);
for j=1:N
    for k=1:S
        Q_S(time_i,j,k)=r.val((time_i-1)*N*S+(j-1)*S+k,4);
    end
end

% CHP energy supply
rs.name = 'E_CHP';
r = rgdx ('results', rs);
E_CHP = zeros(sample_q,N,S);
for j=1:N
    for k=1:S
        E_CHP(time_i,j,k)=r.val((time_i-1)*N*S+(j-1)*S+k,4);
    end
end

% Imbalance reduction
rs.name = 'E_i';
r = rgdx ('results', rs);
E_i = zeros(sample_q,S);
for j=1:S
    E_i(time_i,j)=r.val((time_i-1)*S+j,3);
end

% Bidding price
rs.name = 'E_b';
r = rgdx ('results', rs);
E_b(time_i,1)=r.val(time_i,2);

rs.name = 'ON';
r = rgdx ('results', rs);
ON = zeros(sample_q,N,S);
for j=1:N
    for k=1:S
        ON(time_i,j,k)=r.val((time_i-1)*N*S+(j-1)*S+k,4);
    end
end
end