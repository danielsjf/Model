%% first functional for the Brier scores
% g = 1 for each element of z wich is larger than e on the interval
% [k+h/2,k-h/2]
% g = 0 otherwise

function [g] = functional_val(z,k,h,e)
if h > 1
temp = zeros(h,1);

b = 1; 
if k+h > length(z)
   h = length(z)-k; 
end
for a = k:k+h
   if z(a) >= e
       temp(b) = 1;
   end
   b = b+1;
end
g = prod(temp);
else
    if z(k) >= e
       g = 1; 
    else
        g = 0;
    end
end
end