%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function:     P_versus_A.m
% Calls:        FintgM.m, FintgA.m
% Author:       Mansoor A Haider
% Date:         11/15/2022
% Version:      1.0
% Reference:    MA Haider, KJ Pearce, NC Chesler, NA Hill and MS Olufsen, 
%               Application and reduction of a nonlinear hyperelastic wall 
%               model capturing ex vivo relationships between fluid 
%               pressure, area and wall thickness in normal and hypertensive 
%               murine left pulmonary arteries (under review)
% Purpose:      Function to simulate the pressure ("pressure") as a function 
%               of area (aout) given a fixed set of nominal parameter 
%               values (Rin,cM,cA,k1M,k1A,k2M,k2A,hRef)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function pressure = P_versus_A(aout,Rin,cM,cA,k1M,k1A,k2M,k2A,hRef)

global k Lz hMAFrac mmHgConv

%%% Convert area to radius
rout = sqrt(aout/pi);

%%% Set the media thickness in the stress free state
hMRef=hMAFrac*hRef; 

%%% Simulate the inner radius and the location of the media-adventitia
%%% interface in the current configuration using equation (3)
rin=sqrt(rout^2 - ((Rin+hRef)^2 - Rin^2)/k/Lz);
rMA=sqrt(((Rin+hMRef)^2 - Rin^2)/k/Lz + rin^2);

%%% Via numerical quadtature, approximate the two integrals in equation (9)
pressure = -(integral(@(r)FintgM(r,rin,Rin,cM,k1M,k2M,hRef)*mmHgConv,rin,rMA) + ... 
             integral(@(r)FintgA(r,rin,Rin,cM,cA,k1M,k2A,k1A,hRef)*mmHgConv,rMA,rout));