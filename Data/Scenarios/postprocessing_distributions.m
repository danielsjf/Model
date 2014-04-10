% postprocessing 
clc
clear all
close all

load(['INPUT' filesep 'WPFE_data.mat']); 
load(['INPUT' filesep 'WPFE_ep.mat']);
load(['INPUT' filesep 'WPFE_ee.mat']);

% Prediction errors  --  Forecast
n_int = 1/w_int;
ll_int = [0:w_int:1-w_int];
ul_int = [0+w_int:w_int:1];
ll_error = [-1:w_int:1-w_int];
ul_error = [-1+w_int:w_int:1];
histo_error = zeros(n_int,2*n_int);
histo_prev_error = zeros(2*n_int,2*n_int);

for day = 5:5
load(['SCENARIOS' filesep 'day',num2str(day),'.mat']);
scen_error = scen_red.error;

[samples_scenario,n_scen_gen] = size(scen_error);
forecast = TSwind.forecast((day-1)*samples_scenario+1:day*samples_scenario);
scen_mat = scen_error;
for a = 1:samples_scenario   
    xx = find(forecast(a) >= ll_int & forecast(a) < ul_int);
    if xx >= n_int
        xx = n_int; 
    elseif xx <= 1
        xx = 1; 
    end
    
    for b = 1:n_scen_gen
    yy = find(scen_error(a,b) >= ll_error & scen_error(a,b) < ul_error);
    if yy >= n_int*2
        yy = n_int*2;
    elseif yy <= 1
        yy = 1; 
    end
    
    % histogram (pdf) e-p
    histo_error(xx,yy) = histo_error(xx,yy)+scen_red.prob(b); 
    
    % histogram (pdf) e-e
    if a > 1 
    zz = find(scen_error(a-1,b) >= ll_error & scen_error(a-1,b) < ul_error);
    if zz >= n_int*2
        zz = n_int*2;
    elseif zz <= 1
        zz = 1; 
    end
    histo_prev_error(zz,yy) = histo_prev_error(zz,yy)+scen_red.prob(b);  %+1;%      
    end
    
    end
end

end

% Normalization of the histogram
for a = 1:n_int
   if sum(histo_error(a,:)) ~= 0
   histo_error(a,:) = histo_error(a,:)/sum(histo_error(a,:)); 
   end
end    

for a = 1:2*n_int
   if sum(histo_prev_error(a,:)) ~= 0
   histo_prev_error(a,:) = histo_prev_error(a,:)/sum(histo_prev_error(a,:)); 
   end
end    

histo_cdf_ee = zeros(80,80);
cdf_stable_ee = zeros(80,80);
for a = 1:2*n_int
histo_cdf_ee(a,:) = cumsum(histo_prev_error(a,:));
cdf_stable_ee(a,:) = cumsum(pdf_stable_ee(a,:));
end

histo_cdf_error = zeros(40,80);
cdf_stable = zeros(40,80);
for a = 1:n_int
histo_cdf_error(a,:) = cumsum(histo_error(a,:));
cdf_stable(a,:) = cumsum(pdf_stable(a,:));
end

% ep
error_vec = [-1+w_int/2:w_int:1-w_int/2];
for bin = 1:40
data = [error_vec' pdf_stable(bin,:)' histo_error(bin,:)' cdf_stable(bin,:)' histo_cdf_error(bin,:)'];

if sum(data(:,3)) > 0
figure 
plot(log10(data(:,2:3)));
legend('stable','scenarios')

figure
plot(log10(data(:,4:5)));
legend('stable','scenarios')
end
end

load('WPFE_statistics_ep.mat','histo')
data = [[-1+0.025/2:2/80:1-0.025/2]' [histo.prob(1:36)'*pdf_stable(1:36,:)]' [histo.prob(1:36)'*histo_error(1:36,:)]' [histo.prob(1:36)'*cdf_stable(1:36,:)]' [histo.prob(1:36)'*histo_cdf_error(1:36,:)]'];

figure
plot(data(:,2:3))
legend('overall stable','overall scenarios')

% qq-plots
error = TSwind.production-TSwind.forecast;
y = quantile(error,[0.01:0.01:0.99]);
x = quantile(scen_mat(:),[0.01:0.01:0.99]);
quantiles = [y' x'];

figure
plot(quantiles(:,1),quantiles(:,2),'+')

% ee
% error_vec = [-1+w_int/2:w_int:1-w_int/2];
% 
% for bin = 1:80
% data = [error_vec' pdf_stable_ee(bin,:)' histo_prev_error(bin,:)' cdf_stable_ee(bin,:)' histo_cdf_ee(bin,:)'];
% if sum(data(:,3)) > 0
% 
% figure 
% plot(data(:,2:3));
% legend('stable','scenarios')
% 
% figure
% plot(data(:,4:5));
% legend('stable','scenarios')
% end
% end
% 
% load('WPFE_statistics_ee.mat','histo')
% data2 = [[-1+0.025/2:2/80:1-0.025/2]' [histo.prob'*pdf_stable_ee]' [histo.prob'*histo_prev_error]' [histo.prob'*cdf_stable_ee]' [histo.prob'*histo_cdf_ee]'];
% 
% figure 
% plot(data2(:,2:3));
% legend('stable','scenarios')
% 
% figure
% plot(data2(:,4:5));
% legend('stable','scenarios')
% 
% 
% x2 = [histo.prob'*cdf_stable_ee]';
% x1 = [histo.prob'*histo_cdf_ee]';
% xx = [-1+0.025/2:2/80:1-0.025/2]';
% prob_vec = [0.01:0.01:0.99];
% x_q1 = zeros(length(prob_vec),1);
% x_q2 = zeros(length(prob_vec),1);
% for a= 1:length(prob_vec)
%     x_q1(a) =  xx(find(x1>=prob_vec(a),1,'first'));
%     x_q2(a) =  xx(find(x2>=prob_vec(a),1,'first'));
% end
% 
% quantiles = [x_q2 x_q1];
% 
% figure
% plot(x_q2,x_q1,'+')
% 
