%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function:     FintgM.m
% Author:       Mansoor A Haider
% Date:         11/15/2022
% Version:      1.0
% Reference:    MA Haider, KJ Pearce, NC Chesler, NA Hill and MS Olufsen, 
%               Application and reduction of a nonlinear hyperelastic wall 
%               model capturing ex vivo relationships between fluid 
%               pressure, area and wall thickness in normal and hypertensive 
%               murine left pulmonary arteries (under review)
% Purpose:      Function to generate the integrand for the first integral  
%               in equation (9)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function intg = FintgM(r,rin1,Rin0,cM,k1M,k2M,hRef) 

global k beta Lz

%%% Set the twist angle
Phi=0;
%%% Set the reference axial length at 1 (without loss of generality) when 
%%% the twist angle is zero
l=1;

%%% Use the media fiber angle to calculate the Lagrangian fiber directions 
%%% in equation (6)
a01RM=0; a01HM=cos(beta*pi/180); a01ZM=sin(beta*pi/180);   
a02RM=0; a02HM=cos(beta*pi/180); a02ZM=-sin(beta*pi/180);   

Iv=ones(size(r));
Rin=Rin0*Iv;
rin=rin1*Iv;

%%% Calculate the integrand for the first integral in equation (9)
intgM = 0.1e1./ r.* (cM * (0.2e1./ 0.3e1./ Lz ^ 2./ k ^ 2./ r.^ 2.* (Rin.^ 2 + (r.^ 2 - rin.^ 2) * k * Lz) -...
    r.^ 2.* Phi ^ 2./ l ^ 2 * Lz ^ 2./ 0.3e1 - k ^ 2 * r.^ 2./ (Rin.^ 2 + (r.^ 2 - rin.^ 2) * k * Lz)./ 0.3e1 -...
    Lz ^ 2./ 0.3e1) + 0.2e1 * k1M * (a01RM ^ 2./ Lz ^ 2./ k ^ 2./ r.^ 2.* (Rin.^ 2 + (r.^ 2 - rin.^ 2) * k * Lz) +...
    0.2e1 * a01HM * a01ZM * k * r.^ 2.* (Rin.^ 2 + (r.^ 2 - rin.^ 2) * k * Lz).^ (-0.1e1./ 0.2e1).* Phi./ l * Lz +...
    a01HM ^ 2 * k ^ 2 * r.^ 2./ (Rin.^ 2 + (r.^ 2 - rin.^ 2) * k * Lz) + a01ZM ^ 2 * (r.^ 2.* Phi ^ 2./ l ^ 2 * Lz ^ 2 +...
    Lz ^ 2) - 0.1e1).* exp(k2M * (a01RM ^ 2./ Lz ^ 2./ k ^ 2./ r.^ 2.* (Rin.^ 2 + (r.^ 2 - rin.^ 2) * k * Lz) +...
    0.2e1 * a01HM * a01ZM * k * r.^ 2.* (Rin.^ 2 + (r.^ 2 - rin.^ 2) * k * Lz).^ (-0.1e1./ 0.2e1).* Phi./ l * Lz +...
    a01HM ^ 2 * k ^ 2 * r.^ 2./ (Rin.^ 2 + (r.^ 2 - rin.^ 2) * k * Lz) + a01ZM ^ 2 * (r.^ 2.* Phi ^ 2./ l ^ 2 * Lz ^ 2 +...
    Lz ^ 2) - 0.1e1).^ 2).* (0.2e1./ 0.3e1 * a01RM ^ 2 - a01HM ^ 2./ 0.3e1 - a01ZM ^ 2./ 0.3e1) + 0.2e1 * k1M *...
    (a02RM ^ 2./ Lz ^ 2./ k ^ 2./ r.^ 2.* (Rin.^ 2 + (r.^ 2 - rin.^ 2) * k * Lz) + 0.2e1 * a02HM * a02ZM * k *...
    r.^ 2.* (Rin.^ 2 + (r.^ 2 - rin.^ 2) * k * Lz).^ (-0.1e1./ 0.2e1) * Phi./ l * Lz +...
    a02HM ^ 2 * k ^ 2 * r.^ 2./ (Rin.^ 2 + (r.^ 2 - rin.^ 2) * k * Lz) + a02ZM ^ 2 * (r.^ 2 * Phi ^ 2./ l ^ 2 * Lz ^ 2 +...
    Lz ^ 2) - 0.1e1).* exp(k2M * (a02RM ^ 2./ Lz ^ 2./ k ^ 2./ r.^ 2.* (Rin.^ 2 + (r.^ 2 - rin.^ 2) * k * Lz) +...
    0.2e1 * a02HM * a02ZM * k * r.^ 2.* (Rin.^ 2 + (r.^ 2 - rin.^ 2) * k * Lz).^ (-0.1e1./ 0.2e1).* Phi./ l * Lz +...
    a02HM ^ 2 * k ^ 2 * r.^ 2./ (Rin.^ 2 + (r.^ 2 - rin.^ 2) * k * Lz) + a02ZM ^ 2 * (r.^ 2 * Phi ^ 2./ l ^ 2 * Lz ^ 2 +...
    Lz ^ 2) - 0.1e1).^ 2).* (0.2e1./ 0.3e1 * a02RM ^ 2 - a02HM ^ 2./ 0.3e1 - a02ZM ^ 2./ 0.3e1) - cM * (0.2e1./ 0.3e1 *...
    r.^ 2 * Phi ^ 2./ l ^ 2 * Lz ^ 2 + 0.2e1./ 0.3e1 * k ^ 2 * r.^ 2./ (Rin.^ 2 + (r.^ 2 - rin.^ 2) * k * Lz) -...
    0.1e1./ Lz ^ 2./ k ^ 2./ r.^ 2.* (Rin.^ 2 + (r.^ 2 - rin.^ 2) * k * Lz)./ 0.3e1 - Lz ^ 2./ 0.3e1) - 0.2e1 * k1M *...
    (a01RM ^ 2./ Lz ^ 2./ k ^ 2./ r.^ 2.* (Rin.^ 2 + (r.^ 2 - rin.^ 2) * k * Lz) + 0.2e1 * a01HM * a01ZM *...
    k * r.^ 2.* (Rin.^ 2 + (r.^ 2 - rin.^ 2) * k * Lz).^ (-0.1e1./ 0.2e1) * Phi./ l * Lz + a01HM ^ 2 * k ^ 2.*...
    r.^ 2./ (Rin.^ 2 + (r.^ 2 - rin.^ 2) * k * Lz) + a01ZM ^ 2 * (r.^ 2 * Phi ^ 2./ l ^ 2 * Lz ^ 2 + Lz ^ 2) - 0.1e1).*...
    exp(k2M * (a01RM ^ 2./ Lz ^ 2./ k ^ 2./ r.^ 2.* (Rin.^ 2 + (r.^ 2 - rin.^ 2) * k * Lz) + 0.2e1 * a01HM * a01ZM *...
    k * r.^ 2.* (Rin.^ 2 + (r.^ 2 - rin.^ 2) * k * Lz).^ (-0.1e1./ 0.2e1).* Phi./ l * Lz +...
    a01HM ^ 2 * k ^ 2 * r.^ 2./ (Rin.^ 2 + (r.^ 2 - rin.^ 2) * k * Lz) + a01ZM ^ 2 * (r.^ 2 * Phi ^ 2./ l ^ 2 * Lz ^ 2 +...
    Lz ^ 2) - 0.1e1).^ 2).* (0.2e1./ 0.3e1 * a01HM ^ 2 - a01RM ^ 2./ 0.3e1 - a01ZM ^ 2./ 0.3e1) - 0.2e1 * k1M *...
    (a02RM ^ 2./ Lz ^ 2./ k ^ 2./ r.^ 2.* (Rin.^ 2 + (r.^ 2 - rin.^ 2) * k * Lz) + 0.2e1 * a02HM * a02ZM * k *...
    r.^ 2.* (Rin.^ 2 + (r.^ 2 - rin.^ 2) * k * Lz).^ (-0.1e1./ 0.2e1).* Phi./ l * Lz + a02HM ^ 2 * k ^ 2.*...
    r.^ 2./ (Rin.^ 2 + (r.^ 2 - rin.^ 2) * k * Lz) + a02ZM ^ 2 * (r.^ 2 * Phi ^ 2./ l ^ 2 * Lz ^ 2 + Lz ^ 2) - 0.1e1)...
   .* exp(k2M * (a02RM ^ 2./ Lz ^ 2./ k ^ 2./ r.^ 2.* (Rin.^ 2 + (r.^ 2 - rin.^ 2) * k * Lz) + 0.2e1 * a02HM * a02ZM *...
    k * r.^ 2.* (Rin.^ 2 + (r.^ 2 - rin.^ 2) * k * Lz).^ (-0.1e1./ 0.2e1) * Phi./ l * Lz +...
    a02HM ^ 2 * k ^ 2.* r.^ 2./ (Rin.^ 2 + (r.^ 2 - rin.^ 2) * k * Lz) + a02ZM ^ 2 * (r.^ 2.* Phi ^ 2./ l ^ 2 * Lz ^ 2 +...
    Lz ^ 2) - 0.1e1).^ 2).* (0.2e1./ 0.3e1 * a02HM ^ 2 - a02RM ^ 2./ 0.3e1 - a02ZM ^ 2./ 0.3e1));

intg = intgM;