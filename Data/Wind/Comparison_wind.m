clc;
clear all;

month = 'November';
year = '2013';
date_n = [month,' ',year]; % nice date
date_u = [month,'_',year]; % ugly date

%% Forecasted
file_f = ['Generation_Forecast_Historical_Cipu_',month,'_2013_wind.csv'];
fid = fopen(file_f); %open file
 headers = fgetl(fid);    %get first line
 headers = textscan(headers,'%s','delimiter',';'); %read first line
 format = repmat('%s',1,size(headers{1,1},1)+1); %count columns and make format string
 data = textscan(fid,format,'delimiter',';'); %read rest of the file
 data = [data{:}];
 data_forecast = data(3:end,8:end);
fclose('all');
[j1,j2] = size(data_forecast); j = j1*j2;
data_f = zeros(j1,j2);
tmp1 = zeros(1,j);
i = 1:1:j1;
p = 1:1:j2;
data_forecast = strrep(data_forecast,'(','-');
data_forecast = strrep(data_forecast,')','');
data_f(i,p) = str2double(data_forecast(i,p));
for i=1:j1
    for j=1:j2
        tmp1((i-1)*96+j) = data_f(i,j);
    end
end
data_f = tmp1;

%% Generated
file_g = ['Generation_Produced_Historical_2013-11-',month,'_wind.csv'];
fid = fopen(file_g); %open file
 headers = fgetl(fid);    %get first line
 headers = textscan(headers,'%s','delimiter',';'); %read first line
 format = repmat('%s',1,size(headers{1,1},1)+1); %count columns and make format string
 data = textscan(fid,format,'delimiter',';'); %read rest of the file
 data = [data{:}];
 data_generated = data(4:end,7:end);
fclose('all');
[j1,j2] = size(data_generated); j = j1*j2;
data_g = zeros(j1,j2);
tmp1 = zeros(1,j);
i = 1:1:j1;
p = 1:1:j2;
data_generated = strrep(data_generated,'(','-');
data_generated = strrep(data_generated,')','');
data_g(i,p) = str2double(data_generated(i,p));
for i=1:j1
    for j=1:j2
        tmp1((i-1)*96+j) = data_g(i,j);
    end
end
data_g = tmp1;

%% Imbalance
imbalance = data_f-data_g;
l = length(data_f);
time = (1:l)/96;

figure(1)
plot(time,[imbalance;zeros(1,l)])
title(['Imbalance (',date_n,')']);
xlabel('Time [days]');
ylabel('Imbalance [MW]');
saveas(gcf,['Images\imbalance_',date_u,'.png']);
saveas(gcf,['Images\imbalance_',date_u,'.eps']);

figure(2)
subplot(2,1,1)
plot(time,[imbalance./data_f;zeros(1,l)])
title(['Relative imbalance forcasted (',date_n,')']);
xlabel('Time [days]');
ylabel('Relative imbalance');
subplot(2,1,2)
plot(time,[imbalance./data_g;zeros(1,l)])
title(['Relative imbalance generated (',date_n,')']);
xlabel('Time [days]');
ylabel('Relative imbalance');
saveas(gcf,['Images\relative_imbalance_',date_u,'.png']);
saveas(gcf,['Images\relative_imbalance_',date_u,'.eps']);

%% Imbalance period
tmp1 = [];
tmp2 = [];
tmp3 = 0;
tmp4 = [0 0];
if(imbalance > 0)
    pos = 1;
else
    pos = -1;
end
j = 0;
for i = 1:length(time)
    if((pos*imbalance(i) < 0))
        tmp1 = [tmp1 i-j];
        tmp2 = [tmp2 tmp3/(i-j)];
        p = 1.5 + pos*0.5;
        tmp4(p) = tmp4(p) + tmp3;
        tmp3 = 0;
        pos = -pos;
        j = i;
    end
    tmp3 = tmp3 + imbalance(i);
end
imbal_cumul_hours = tmp1/4;
imbal_avg_power = tmp2;

figure(3)
plot(imbal_cumul_hours)
title(['Periods of under -or overproduction (',date_n,')']);
xlabel('Periods');
ylabel('Length [hours]');

spread = 1/4:1/4:max(imbal_cumul_hours);
spread2 = zeros(size(spread));
for i = spread*4
    tmp = find(imbal_cumul_hours==spread(i));
    spread2(i) = length(tmp);
end

figure(4)
plot(spread,spread2)
title(['Periods of under -or overproduction (',date_n,')']);
xlabel('Length [hours]');
ylabel('Frequency');

figure(5)
plot(imbal_avg_power)
title(['Periods of under -or overproduction (',date_n,')']);
xlabel('Periods');
ylabel('Average imbalance [MW]');

powerstep = 5;
tmp2 = abs(imbal_avg_power);
spread = 0:powerstep:ceil(max(tmp2)/powerstep)*powerstep;
spread2 = zeros(size(spread));
for i = 1:length(spread(1:end-1));
    tmp = find(tmp2(find(tmp2>spread(i)))<=spread(i+1));
    spread2(i) = length(tmp);
end

figure(6)
plot(spread,spread2)
title(['Periods of under -or overproduction (',date_n,')']);
xlabel('Average imbalance [MW]');
ylabel('Frequency');

span = 2*96;
figure(7)
subplot(2,1,1)
plot(time(1:span),[data_f(1:span);data_g(1:span)])
title('Forecast VS generated');
xlabel('Time [days]');
ylabel('Power [MW]');
legend('Forecast','Generated');

subplot(2,1,2)
plot(time(1:span),[imbalance(1:span);zeros(1,span)])
title('Imbalance (=error)');
xlabel('Time [days]');
ylabel('Power [MW]');
saveas(gcf,['Images\ForecastVSGenerated.png']);
saveas(gcf,['Images\ForecastVSGenerated.eps']);

bidding = 150;
figure(8)
plot(time(1:span),[max(ones(1,span)*bidding-imbalance(1:span),0);ones(1,span)*bidding])
title('Bidding');
xlabel('Time [days]');
ylabel('Power [MW]');
ylabel('Power [MW]');
legend('Bidding minus error','Bidding');
saveas(gcf,['Images\Bidding.png']);
saveas(gcf,['Images\Bidding.eps']);

span = 0.5*96;
figure(9)
plot(time(1:span),[imbalance(1:span);imbalance(1:span)+(rand(1,span)-0.5)*4*mean(abs(imbalance(1:span)));(1:span)*0.1+(rand(1,span)-0.5)*4*mean(abs(imbalance(1:span)))]);
title('Probability');
xlabel('Time [days]');
ylabel('Power [MW]');
ylabel('Power [MW]');
legend('Scenario 1: P=0.25','Scenario 2: P=0.30','Scenario 3: P=0.45');
saveas(gcf,['Images\BiddingProb.png']);
saveas(gcf,['Images\BiddingProb.eps']);

span = 2*96;
figure(10)
plot(time(1:span),[data_f(1:span);data_g(1:span)])
title('Forecast VS generated');
xlabel('Time [days]');
ylabel('Power [MW]');
legend('Forecast','Generated');
saveas(gcf,['Images\ForecastVSGenerated2.png']);
saveas(gcf,['Images\ForecastVSGenerated2.eps']);