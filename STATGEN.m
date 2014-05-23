function [GEN] = STATGEN(DATAsample,DATAgen,dib,Period)
if nargin < 3
  dib = 0;
end
if nargin < 4
  Period = 1:size(DATAsample.heatD_s,1);
end

GEN.costWOchpT = 0;
for i=1:size(DATAsample.heatD_s,2)
    GEN.costWOchpT = GEN.costWOchpT - sum(DATAsample.heatD_s(Period,i) .* DATAsample.price_gas_s(Period) / DATAgen.Nb0); % Cost without CHP for one unit
end

if dib == 1, disp(['Cost without CHP: ',num2str(GEN.costWOchpT)]); end
end