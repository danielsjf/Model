clc
close all

set(0,'DefaultAxesFontName', 'Times')
set(0,'DefaultAxesFontSize', 9)
set(0,'DefaultUiControlUnits','normalized');
set(0,'DefaultFigurePaperUnits', 'centimeters');
set(0,'DefaultFigurePaperPosition', [10 10 14 10]);
LFS = 7; % Legend font size

% Run CHP.m before running this script

CHPunits % CHP data

temp = pwd;
ThesisRoot = temp(1:numel(temp)-17);
ThesisImages = [ThesisRoot, filesep, 'Thesis', filesep, 'Images', filesep];

load Load2013.mat

figure()
plot(time_q,[elecWF_s,elecWP_s])
title('Wind production forecast and actual production');
legend({'Wind production forecast','Wind production'},'FontSize',LFS);
xlabel('Time [h]');
ylabel('Energy [MWh]');
saveas(gcf,[ThesisImages,'Wind_ForecastAndProduction.png']);
print('-depsc2',[ThesisImages,'Wind_ForecastAndProduction.eps']);

figure()
subplot(2,1,1)
%set(gcf,'DefaultAxesColorOrder',[0.8 0 0.8; 0.8 0.6 0])
plot(time_q,actuals_s-0.05, '-', 'Color', [0.8 0 0.8])
hold on
plot(time_q,zeros(size(actuals_s)), '-', 'Color', [0.8 0.6 0])
axis([1 48 -0.1 0.2]);
title({'Wind forecast error'});
xlabel('Time [h]');
ylabel('Forecast error [MWh]');
subplot(2,1,2)
bidIC = 0.8*max(actuals_s-0.05)*ones(size(actuals_s)); % Installed capacity
bidMC = 0.3*bidIC; % Minimal capacity
bidRB = 0.6*bidIC; % Real bid
bidP = max(min(bidRB-actuals_s+0.05,bidIC),bidMC);
bidP(bidP == bidMC) = 0;
reset(gca)
plot(time_q,bidIC, '--', 'Color', [0.8 0 0])
hold on
plot(time_q,bidMC, '--', 'Color', [0 0.6 0])
hold on
plot(time_q,bidRB, '-.', 'Color', [0 0 0.8])
hold on
plot(time_q,bidRB-actuals_s+0.05, '-', 'Color', [0.8 0 0.8])
hold on
plot(time_q,bidP, '-', 'Color', [0 0.6 0.8])
axis([1 48 -0.1 0.2]);
title('CHP bid and production in VPP');
legend({'CHP inst. cap.','CHP min. cap.','CHP bid','CHP opt. prod.','CHP real prod.'},'Location','Southwest','FontSize',LFS);
xlabel('Time [h]');
ylabel('Energy [MWh]');
saveas(gcf,[ThesisImages,'CHP_BidAndProduction.png']);
print('-depsc2',[ThesisImages,'CHP_BidAndProduction.eps']);

figure()
plot(time_q,actuals_s)
hold on
area(time_q,-actuals_s, 'FaceColor', [0 0.8 0], 'EdgeColor', [0 0.3 0])
axis([1 48 -0.2 0.2]);
title('Wind forecast error and imbalance reduction zone');
legend({'Wind forecast error (imbalance)','Imbalance reduction zone'},'FontSize',LFS);
xlabel('Time [h]');
ylabel('Energy [MWh]');
saveas(gcf,[ThesisImages,'Wind_ImbalanceAndReduction.png']);
print('-depsc2',[ThesisImages,'Wind_ImbalanceAndReduction.eps']);

figure()
plot(price_elecS_y,POS,'b.',price_elecS_y,NEG,'r.')
title('Relation spot price - imbalance tariff');
legend({'POS','NEG'},'FontSize',LFS);
xlabel('Spot price [€/MWh]');
ylabel('Imbalance tariff [€/MWh]');
saveas(gcf,[ThesisImages,'ImbalP_POSNEGSpot.png']);
print('-depsc2',[ThesisImages,'ImbalP_POSNEGSpot.eps']);

figure()
Day = reshape(price_elecS_y, 96, 365)';
plotshaded(1:96,Day,'r',1,0.1)
title('Spot price (extrema and mean per quarter)');
xlabel('Time [quarters]');
ylabel('Spot price [€/MWh]');
saveas(gcf,[ThesisImages,'Spot_Yearly.png']);
print('-depsc2',[ThesisImages,'Spot_Yearly.eps']);

figure()
Day = reshape(price_elecS_y, 96, 365)';
plotshaded(1:365,Day','r',1,0.1)
title('Spot price (extrema and mean per day)');
xlabel('Time [days]');
ylabel('Spot price [€/MWh]');
saveas(gcf,[ThesisImages,'Spot_Daily.png']);
print('-depsc2',[ThesisImages,'Spot_Daily.eps']);

figure()
subplot(2,1,1)
Day = reshape(POS, 96, 365)';
plotshaded(1:365,Day','b',1,0.1)
Day = reshape(price_elecS_y, 96, 365)';
plotshaded(1:365,Day','r',1,0.1)
title('Spot price and imbalance tariff (POS)');
xlabel('Time [days]');
ylabel('Spot price [€/MWh]')
subplot(2,1,2)
Day = reshape(NEG, 96, 365)';
plotshaded(1:365,Day','b',1,0.1)
Day = reshape(price_elecS_y, 96, 365)';
plotshaded(1:365,Day','r',1,0.1)
title('Spot price and imbalance tariff (NEG)');
xlabel('Time [days]');
ylabel('Spot price [€/MWh]')
saveas(gcf,[ThesisImages,'Spot_POSNEG.png']);
print('-depsc2',[ThesisImages,'Spot_POSNEG.eps']);

figure()
plot(NRV,POS - price_elecS_y,'b.',NRV,NEG - price_elecS_y,'r.')
title('Relation NRV - (imbalance tariff - spot price)');
legend({'POS - spot price','NEG - spot price'},'FontSize',LFS);
xlabel('NRV [MWh]');
ylabel('Imbalance tariff - spot price [€/MWh]');
saveas(gcf,[ThesisImages,'ImbalP_POSNEGSpotNRV.png']);
print('-depsc2',[ThesisImages,'ImbalP_POSNEGSpotNRV.eps']);

figure()
plot(NRV,POS,'b.',NRV,NEG,'r.')
title('Relation NRV - imbalance tariff');
legend({'POS','NEG'},'FontSize',LFS);
xlabel('NRV [MWh]');
ylabel('Imbalance tariff [€/MWh]');
saveas(gcf,[ThesisImages,'ImbalP_POSNEGNRV.png']);
print('-depsc2',[ThesisImages,'ImbalP_POSNEGNRV.eps']);

figure()
plot([price_elecS_y,POS,NEG,NRV]) %,NRV,NEG,'r.')
title('NRV, spot price, imbalance tariff');
legend({'Spot price','POS','NEG','NRV'},'FontSize',LFS);
xlabel('Time [quarters]');
ylabel('Price [€/MWh]');
saveas(gcf,[ThesisImages,'ImbalP_POS-NEG-Spot-NRV.png']);
print('-depsc2',[ThesisImages,'ImbalP_POS-NEG-Spot-NRV.eps']);

figure()
plot([POS-price_elecS_y,NRV]) %,NRV,NEG,'r.')
title('Imbalance tariff - spot price, NRV');
legend({'Imbalance tariff - spot price','NRV'},'FontSize',LFS);
xlabel('Time [quarters]');
ylabel('Price [€/MWh]');
saveas(gcf,[ThesisImages,'ImbalP_POSSpot-NRV.png']);
print('-depsc2',[ThesisImages,'ImbalP_POSSpot-NRV.eps']);

figure()
NRVMDP = fit(NRV(NRV<0),MDP(NRV<0)-price_elecS_y(NRV<0),'poly1')
plot(NRVMDP,'k-',NRV(NRV<0),MDP(NRV<0)-price_elecS_y(NRV<0),'r.')
hold on
NRVMIP = fit(NRV(NRV>0),MIP(NRV>0)-price_elecS_y(NRV>0),'poly1')
plot(NRVMIP,'k--',NRV(NRV>0),MIP(NRV>0)-price_elecS_y(NRV>0),'b.')
title('MDP and MIP')
xlabel('NRV [MWh]');
ylabel('Imbalance tariff [€/MWh]');
legend({'MDP - spot price','MDP - spot price (fit)','MIP - spot price','MIP - spot price (fit)'},'FontSize',LFS)
saveas(gcf,[ThesisImages,'ImbalP_MDPMIPfit.png']);
print('-depsc2',[ThesisImages,'ImbalP_MDPMIPfit.eps']);

figure()
plot(abs(NRV(NRV>0)),alpha(NRV>0),'.','Color',[1,0.5,0.5])
hold on
plot(abs(NRV(NRV<0)),alpha(NRV<0),'.','Color',[0.5,0.5,1])
NRVAlphaPos = fit(abs(NRV(NRV>0)),alpha(NRV>0),'poly1')
plot(NRVAlphaPos,'r-.')
NRVAlphaNeg = fit(abs(NRV(NRV<0)),alpha(NRV<0),'poly1')
plot(NRVAlphaNeg,'b--')
NRVAlpha = fit(abs(NRV),alpha,'poly1')    
plot(NRVAlpha,'k-')
axis([0 max(NRV) 0 max(alpha)])
legend({'Alpha (NRV > 0)','Alpha (NRV < 0)','Alpha (NRV > 0) (fit)','Alpha (NRV > 0) (fit)','Alpha (fit)'},'FontSize',LFS)
title('Alpha as a function of NRV')
xlabel('|NRV| [MWh]');
ylabel('Alpha [€/MWh]');
saveas(gcf,[ThesisImages,'Alpha.png']);
print('-depsc2',[ThesisImages,'Alpha.eps']);

figure()
for i = 1:3
    subplot(3,1,i)
    plot(1/96:1/96:365,heatD_y(:,i)*1000, 'b', 1/96:1/96:365, mean(heatD_y(:,i)*1000)*ones(size(heatD_y,1),1),'r')
    axis([0 365 0 max(max(heatD_y))*1000])
    title(['Heat demand for unit ',num2str(i)])
    xlabel('Time [days]');
    ylabel({'Heat demand';'[kWh]'});
end
saveas(gcf,[ThesisImages,'HeatDemand1.png']);
print('-depsc2',[ThesisImages,'HeatDemand1.eps']);

figure()
for i = 4:6
    subplot(3,1,i-3)
    plot(1/96:1/96:365,heatD_y(:,i)*1000, 'b', 1/96:1/96:365, mean(heatD_y(:,i)*1000)*ones(size(heatD_y,1),1),'r')
    axis([0 365 0 max(max(heatD_y))*1000])
    title(['Heat demand for unit ',num2str(i)])
    xlabel('Time [days]');
    ylabel({'Heat demand';'[kWh]'});
end
saveas(gcf,[ThesisImages,'HeatDemand2.png']);
print('-depsc2',[ThesisImages,'HeatDemand2.eps']);

figure()
hold on
p1=plot(duration(:,1),thermalload(:,1));
h2=rectangle('Position',[0,0,duration_opt(1),thermalload_opt(1)]);
p2=plot(nan,nan,'s','markeredgecolor',get(h2,'edgecolor'));
legend([p1,p2],{'Thermal load','Largest rectangle'},'FontSize',LFS);
title(['Heat-load duration diagram for unit ',num2str(1)]);
xlabel('Load duration [h]');
ylabel('Thermal load [MW_t]');
saveas(gcf,[ThesisImages,'House_LoadDurDia.png']);
print('-depsc2',[ThesisImages,'House_LoadDurDia.eps']);

figure()
plot(1/96:1/96:29,heatD_y(31*96+1:60*96,1)*1000, 'b', 1/96:1/96:29,elecD_y(31*96+1:60*96,1)*1000, 'r')
axis([1/96 29 0 max(max(heatD_y(31*96+1:60*96,1)))*1000*1.1])
title('Heat and electricity demand')
legend({'Electricity demand','Thermal demand'},'FontSize',LFS)
xlabel('Time [days]');
ylabel({'Heat demand';'[kWh]'});
saveas(gcf,[ThesisImages,'House_QandPdemand.png']);
print('-depsc2',[ThesisImages,'House_QandPdemand.eps']);

cogen(:,[2,4])
color = ['r.'; 'g.'; 'k.'; 'm.'; 'b.'; 'c.'];
figure()
for i = 1:6
    temp = cogen(:,1);
    temp = find(temp==i);
    cogen(temp,:)
    plot(cogen(temp,2),cogen(temp,4),color(i,:))
    hold on
end
legend({'ICE','Stirling','Fuel cell (SOFC)','Rankine','Fuel cell (PEM)'},'FontSize',LFS)
saveas(gcf,[ThesisImages,'CHP_unitdata.png']);
print('-depsc2',[ThesisImages,'CHP_unitdata.eps']);