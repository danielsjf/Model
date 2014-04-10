time_period = 0.25;
nnodes = 10;

gen = xlsread('DATA_GENERATION.xlsx','conventional park','D3:Q73');
ngen = length(gen(:,1));
minup = gen(:,13)/time_period;
mindown = gen(:,14)/time_period;
gen_max = gen(:,2); 
gen_int = gen(:,3);
gen_min = gen(:,4); 
startup_cost = gen(:,6); 
delta_max_up = gen(:,11)*time_period*60;
delta_max_down = gen(:,12)*time_period*60;
for a = 1:ngen
   if delta_max_up(a) > gen_max(a)-gen_min(a)
      delta_max_up(a) = gen_max(a)-gen_min(a);
      delta_max_down(a) = gen_max(a)-gen_min(a);
   end
end
% location characteristics
location = zeros(ngen,nnodes);
location(:,gen(:,1)) = 1;
% fuel prices
fuel_price = gen(:,5);
% Partial loading characteristics
C = fuel_price./gen(:,10).*gen_min;
MA = (fuel_price./gen(:,8).*gen_max-fuel_price./gen(:,10).*gen_min)./(gen_max-gen_min);
D = gen(:,7).*gen_min./(gen(:,10)./gen(:,8));
MB = (gen(:,7).*gen_max-gen(:,7)./(gen(:,10)./gen(:,8)).*gen_min)./(gen_max-gen_min);
