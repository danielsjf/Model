%% first functional for the Brier scores
% g = 1 for each element of z wich is larger than e on the interval
% [k+h/2,k-h/2]
% g = 0 otherwise

function [g] = functional_var(z,k,h,e)
if k+h > length(z)
   h = length(z)-k; 
end

temp_max = max(z(k:k+h));
temp_min = min(z(k:k+h));
diff = temp_max-temp_min;

if diff >= e
   g = 1;
else 
   g = 0;
end
end