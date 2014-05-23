function [Bool] = TEST(DATAsample,DATAres,dib,err)
if nargin < 4
  err = 0;
end
if nargin < 3
  dib = 0;
end

% Heat
EqHeat = sum(DATAres.Q_CHP(:,:,1),2) + sum(DATAres.DeltaQ_S(:,:,1),2) + sum(DATAres.Q_B(:,:,1),2) - sum(DATAsample.heatD_s(:,:,1),2);
if dib == 1; figure(); plot(EqHeat); end
if err == 1 && max(abs(EqHeat)) > 1e-15; error('infeasible'); end

% Power
EqPower = sum(DATAres.Q_CHP(:,:,1),2) + sum(DATAres.DeltaQ_S(:,:,1),2) + sum(DATAres.Q_B(:,:,1),2) - sum(DATAsample.heatD_s(:,:,1),2);
if dib == 1; figure(); plot(EqHeat); end
if err == 1 && max(abs(EqHeat)) > 1e-15; error('infeasible'); end
end