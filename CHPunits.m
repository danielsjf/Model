close all;
clc;

% First column: type (1 = ICE; 2 = Stirling; 3 = Turbine; 4 = Fuel cell (SOFC); 5 = Turbine (Rankine cycle); 6 = Fuel cell (PEM); 7 = Turbine (Brayton cycle))
% Second column: condensor (1 = Y; 2 = N)
% Second column: kWe
% Third column: kWt
% Fourth column: alphaE (LHV)
% Fifth column: alphaQ (LHV)
cogen = [
    2, 2, 1, 7, 0.11, 0.77; % WhisperGen
    2, 2, 1, 3, 0.20, 0.60; % MicroGen
    2, 2, 1, 4, 0.17, 0.68; % Infinia
    2, 2, 3, 15, 0.16, 0.80; % Disenco (inspirit)
    2, 2, 9, 26, 0.25, 0.72; % Cleanergy
    2, 2, 5, 19, 0.20, 0.75; % Qnergy
    5, 2, 1, 10, 0.07, 0.70; % Genlec
    5, 2, 3, 30, 0.07, 0.70; % Climate energy
    5, 2, 2, 16, 0.08, 0.68; % Otag
    5, 2, 2.5, 12, 0.13, 0.65; % Cogen Micro
    1, 2, 1.2, 3, 0.24, 0.60; % Ecowill
    1, 2, 1, 2.5, 0.263, 0.657; % Vaillant EcoPower 1.0
    1, 2, 2, 5.5, 0.25, 0.69; % Proenvis
    1, 2, 1.9, 9, 0.19, 0.76; % Kirsch nano
    1, 2, 5.5, 13.5, 0.27, 0.67; % Viessmann
    1, 2, 6, 13.5, 0.265, 0.595; % Aisin Seiki
    1, 2, 5.5, 15.5, 0.269, 0.611; % Baxi Dachs Senertec
    1, 1, 5.5, 15.5, 0.269, 0.758; % Baxi Dachs Senertec (with condensor)
    1, 2, 3.0, 8.0, 0.25, 0.665; % Vaillant EcoPower 3.0
    1, 2, 4.7, 12.5, 0.25, 0.665; % Vaillant EcoPower 4.7
    1, 2, 2, 8, 0.25, 0.70; % Kirsch micro
    1, 2, 3, 10, 0.25, 0.70; % Kirsch micro
    1, 2, 4, 12, 0.25, 0.70; % Kirsch micro
    1, 2, 3.8, 10.7, 0.24, 0.68; % Proenvis
    1, 2, 3.87, 8.38, 0.267, 0.578; % Yanmar CP4WE
    1, 2, 5.0, 13.0, 0.263, 0.657; %RMB energie
    1, 2, 6.0, 13.5, 0.295, 0.635; % EC Power XRGI 6 (without condenser) (ook dB info verkrijgbaar)
    1, 2, 9.0, 20.0, 0.295, 0.635; % EC Power XRGI 9 (without condenser) (ook dB info verkrijgbaar)
    1, 2, 15.0, 30.0, 0.30, 0.62; % EC Power XRGI 15 (without condenser) (ook dB info verkrijgbaar)
    1, 2, 20, 40, 0.32, 0.64; % EC Power XRGI 20 (without condenser) (ook dB info verkrijgbaar)
    1, 2, 10.0, 17.4, 0.307, 0.533; % Yanmar CP10WE
    1, 2, 25.0, 38.4, 0.335, 0.515; % Yanmar CP25WE-TN
    1, 2, 19.0, 36.0, 0.329, 0.618; % LichtBlick (without condenser)
    1, 1, 19.0, 41.0, 0.329, 0.704; % LichtBlick (with condenser)
    1, 2, 20.0, 42.0, 0.295, 0.620; % Valliant EcoPower 20.0
    1, 2, 30.0, 53.2, 0.320, 0.568; % Kirsch mini30
    4, 2, 1, 1.8, 0.339, 0.611; % Hexis Galileo
    4, 2, 1.5, 0.6, 0.60, 0.250; % Ceramic Fuel Cells BlueGen
    4, 2, 2, 1, 0.57, 0.285; % Ceramic Fuel Cells Gennex
    4, 2, 1, 1.7, 0.30, 0.51; % Vaillant
    4, 2, 0.7, 0.65, 0.465, 0.435; % Kyocera ENE-farm type S
    4, 2, 1, 0.9, 0.45, 0.85; % Topsoe
    4, 2, 1, 2, 0.31, 0.6; % SOFC power Engen 500
    6, 2, 0.3, 0.6, 0.327, 0.653; % Elcore 2400
    6, 2, 5, 7.5, 0.34, 0.58; % Inhouse-engineering Inhouse5000+
    6, 2, 0.75, 1, 0.37, 0.53; % Panasonic/Viessmann (4000 cycli)
    7, 2, 3.0, 14.4, 0.15, 0.72] % Enertwin
    7, 2, 28, 0.25, 
;

cogen = [cogen sum(cogen(:,5:6),2)];

xlswrite('CHPunits.xlsx',cogen);

for i=1:7
    temp1 = find(cogen(:,1)==i);
    temp1 = cogen(temp1,:);
    temp2 = find(temp1(:,2)==2); % no condensor
    temp3 = find(temp1(:,2)==1); % with condensor
    disp('type     minE     maxE     minT     maxT');
    if numel(temp2) > 0; disp([num2str(i), '        ', num2str(min(temp1(temp2,5)),'%1.2f'), '     ', num2str(max(temp1(temp2,5)),'%1.2f'), '     ', num2str(min(temp1(temp2,7)),'%1.2f'), '     ', num2str(max(temp1(temp2,7)),'%1.2f')]); end
    if numel(temp3) > 0; disp([num2str(i), '        ', num2str(min(temp1(temp3,5)),'%1.2f'), '     ', num2str(max(temp1(temp3,5)),'%1.2f'), '     ', num2str(min(temp1(temp3,7)),'%1.2f'), '     ', num2str(max(temp1(temp3,7)),'%1.2f')]); end
    disp(' ');
end

% Source: http://www.microchap.info/ and company websites