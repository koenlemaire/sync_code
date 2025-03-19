function [output]=inv_jumper(t,state,phidd,basedd,grf,grm,parms)
%function [statedot,output]=bouncingstickbasefixed(t,state,u0,parms)
%
% deze functie kan worden aangeroepen door een standaard integrator;
% gebruikmakend van segdyn (Casius et al., 2004) wordt de versnelling van 
% een keten van segmenten berekend
%
% inputs: 
% t (scalar, huidige tijd)
% state (vector, huidige waarde van de toestand van het systeem; zie paragraaf 5.5 van Casius et al.)
% u0 (vector, steady state waarden voor de ingangen van het systeem)
% parms (struct, bevat alle verder benodigde parameters, enige eis is dat
% de velden binnen deze struct in de verschillende functies consistent met
% elkaar zijn
%
% outputs:
% statedot (KOLOMvector, afgeleide van state en dus zelfde lengte als state)
% y (KOLOMvector, alle uitgangen van het systeem)
%   NB y wordt genegeerd tijdens het uitvoeren van een simulatie met een standaard MATLAB integrator
%   daarom is het ivm rekentijd verstandig y alleen te berekenen als dat
%   nodig is; dit is te regelen vanuit het hoofdprogramma via
%   parms.calculate_outputs; zie voorbeeld hoofdprogramma

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% applicatiespecifiek: segdynstate (volledige vector coordinaten (zie segdyn)) maken uit state
% in veel gevallen zullen segdynstate en state identiek zijn
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

segdynstate=state;
L=parms.segparms.L;
d=parms.segparms.d;
m=parms.segparms.m;
j=parms.segparms.j;
nseg=length(L);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% applicatiespecifiek: uitwendige krachten en momenten (voor zover niet onbekend) berekenen
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Fextx=zeros(1,nseg);
Fexty=-9.81*m(:)';
Mext=zeros(1,nseg);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% applicatiespecifiek: (optioneel) tijds- of toestands-afhankelijk deel van u berekenen
% bij constante input zal u gelijk zijn aan u0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% applicatiespecifiek: K en V aanmaken (voor beschrijving: zie segdyn)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

K=...
    [1 0 0 0 0 ...% Frx (n+1)
     1 0 0 0 0 ...% Fry (n+1)
     1 0 0 0 0 ...% M (n+1)
     1 1 1 1 ...% Fextx (n)
     1 1 1 1 ...% Fexty (n)
     1 1 1 1 ...% Mext (n)     
     1 1 1 1 ...% fidd (n)
     1 1 ];% basedd (2)
V=...
    [grf(1) nan nan nan nan ...% Frx (n+1)
     grf(2) nan nan nan nan ...% Fry (n+1)
     grm nan nan nan nan ...% M (n+1)
     Fextx ...% Fextx (n)
     Fexty ...% Fexty (n)
     Mext  ...% Mext (n)     
     phidd ...% fidd (n)
     basedd ];% basedd (2)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% segdyn aanroepen
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[segdynstatedot, Vnew, succes] = segdyn(segdynstate, parms.segparms, K, V);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% applicatiespecifiek: statedot maken uit segdynstatedot
% in veel gevallen zullen segdynstatedot en statedot identiek zijn
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% applicatiespecifiek: uitgangen y berekenen
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% y is een KOLOMvector van willekeurige lengte die alle interessant geachte afhankelijke variabelen
% bevat (bijv reactiekrachten, versnellingen, energietermen, zwaartepuntpositie etc etc)

output=Vnew(:);
