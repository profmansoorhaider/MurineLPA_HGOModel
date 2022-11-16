%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script:       HGOWallExample.m
% Calls:        P_versus_A.m
% Author:       Mansoor A Haider
% Date:         11/15/2022
% Version:      1.0
% Reference:    MA Haider, KJ Pearce, NC Chesler, NA Hill and MS Olufsen, 
%               Application and reduction of a nonlinear hyperelastic wall 
%               model capturing ex vivo relationships between fluid 
%               pressure, area and wall thickness in normal and hypertensive 
%               murine left pulmonary arteries (under review)
% Purpose:      Main script to run the forward model and predict relations
%               between pressure and area, and wall thickness and area for 
%               a fixed set of nominal parameter values.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all

%%% Declare fixed global parameters needed in functions called by this script

global beta     % collagen fiber angle in media and adventitia (Table 4)
global Lz       % axial stretch in deformation (Table 4)
global k        % dimensionless parameter involving opening angle (Eq. 3)
global hMRef    % vessel wall thickness in stress-free state (Table 4) 
global hMAFrac  % fraction of vessel wall thickness occupied by the media
global mmHgConv % conversion factor for pressure units

mmHgConv=0.75006375541921e-2; 

%%% Pressure and outer diameter data from Tabima et al 2010 (Ref [27])
PData2010 = [0 5 10 15 20 25 30 35 40 45 50]';
ODData2010_CTL = [530.4 581.1 662.8 791.1 966.0 1080.3  1130.1 1153.7 1167.6 1178.1 1186.7]';
ODData2010_HPH = [523.2 584.3 637.4 691.7 747.1 790.5 821.9 844.7 863.7 877.3 889.5]';

%%% Wall thickness data from Tabima et al 2010 (Ref [27])
WT_CTL_Mean = 1e-6*[31.33 16.18 15.70]; % for p = 10, 30 & 40 mmHg
WT_HPH_Mean = 1e-6*[45.90 35.16 34.13]; % for p = 10, 30 & 40 mmHg

%%% Choose the data sets for the healthy (CTL) case
WT_DATA = WT_CTL_Mean;
ODData2010 = 1e-6*ODData2010_CTL;

%%% Set the collagen fiber angle based on the reported value in Ref [33]
%%% (see Table 4)
beta = 90 - 35.55;

%%% Set the fraction of vessel wall thickness occupied by the media based 
%%% on the reported value in Ref [33](see Table 4). 
hMAFrac = 0.63;

%%% Set the opening angle based on the reported values in Refs. [25,26]
alpha=94.2*pi/180;

%%% Dimensionless parameter involving opening angle (Eq. 3)
k=2*pi/(2*pi-alpha);

%%% Set the axial stretch from Tabima et al 2010 (Ref [27]), see also 
%%% Sec. 2.5
Lz = 1.4;

%%% Convert outer diameter data to data for outer area and outer radius
NData=11;
AMouse = pi*(ODData2010(1:NData)./2).^2;
routData=sqrt(AMouse/pi);

%%% Determine minimum and maximum outer radius in the data
routMin=min(routData);
routMax=max(routData);

%%% Nominal parameter values for Baseline case with m=8 (Table 1)
RinNom=3.767105750989544e-04;   % in m
cMNom=2.483472503217928e+04;    % in Pa
k1MNom=2.704723554837858e+02;   % in Pa
k2MNom=2.064299307431537;       % dimensionless   
k2ANom=2.734100045646786;       % dimensionless
cANom=1.646039803327283e+04;    % in Pa
k1ANom=34.591557405547093;      % in Pa
hRefNom=4.542755757382240e-05;  % in m

%%% Set the media thickness in the stress free state
hMRef=hMAFrac*hRefNom; 

%%% Initialize vectors to store model predictions of inner radius, 
%%% media-adventitia interface location and wall thickness in the 
%%% current configuration
rinModel=zeros(NData,1);
rMAModel=zeros(NData,1);
WTModel=zeros(NData,1);

%%% Simulate the inner radius and the location of the media-adventitia
%%% interface in the current configuration at the data points using 
%%% equation (3)
for ii=1:NData
    rinModel(ii)=sqrt(routData(ii)^2 - ((RinNom+hRefNom)^2 - RinNom^2)/k/Lz);
    rMAModel(ii)=sqrt(((RinNom+hMRef)^2 - RinNom^2)/k/Lz + rinModel(ii)^2);
    WTModel(ii)=routData(ii) - rinModel(ii);
end

%%% Set resolution and initialize vector for plotting pressure-area
%%% relations
NPlot=100;
pPlot=zeros(NPlot,1);

%%% Generate simulated relations between pressure and the inner area and
%%% the outer area
for i=1:NPlot
    routVec(i)=routMin + (i-1)*(routMax - routMin)/(NPlot-1);
    aout(i)=pi*routVec(i)^2;
    %%% Simulate inner radius via equation (3)
    rinSim(i)=sqrt(routVec(i)^2 - ((RinNom+hRefNom)^2 - RinNom^2)/k/Lz);
    ainSim(i)=pi*rinSim(i)^2;
    %%% Predict the pressure using the pressure-area relation (9)
    pPlot(i) = P_versus_A(aout(i),RinNom,cMNom,cANom,k1MNom,k1ANom,k2MNom,k2ANom,hRefNom);
end

%%% Plot the pressure-area relations and the data
figure(1)
scaleFac=1e12;
plot(aout*scaleFac,pPlot,'k','LineWidth',3);
hold on
plot(AMouse*scaleFac,PData2010,'ko','MarkerSize',14);
plot(ainSim*scaleFac,pPlot,'r','LineWidth',3);
xlabel('Area (um^2)')
ylabel('Pressure (mmHg)')
legend('Outer Area','Data','Inner Area')
xlim([1e-7*scaleFac 12e-7*scaleFac])
ylim([0 55])
ax = gca;
ax.FontSize = 20;

%%% Plot the wall thickness vs pressure relation and the data
figure(2)
plot(PData2010,WTModel*1e6,'ko','LineWidth',3);
hold on
plot([10 30 40],WT_DATA*1e6,'ks','MarkerSize',14);
xlabel('Pressure (mmHg)','FontSize',20)
ylabel('Wall Thickness (um)','FontSize',20)
set(gca,'FontSize',18)
legend('Model','Data')