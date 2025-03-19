function [statedot,output]=frwd_jumper(t,state,parms)
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


L=parms.segparms.L;
d=parms.segparms.d;
m=parms.segparms.m;
j=parms.segparms.j;
nseg=length(L);

state=state(:); % make column vector
segdynstate=state;
phi=state(1:nseg);
phid=state(nseg+1:2*nseg);
base=state(2*nseg+1:2*nseg+2);
based=state(2*nseg+3:2*nseg+4);

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

control=parms.controls;
t_control=parms.t_control;
Mjo = interp1(t_control,control',t);

%Mjo=inputTorque(t,state);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% applicatiespecifiek: K en V aanmaken (voor beschrijving: zie segdyn)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

K=...
    [0 0 0 0 1 ...% Frx (n+1)
     0 0 0 0 1 ...% Fry (n+1)
     1 1 1 1 1 ...% M (n+1)
     1 1 1 1 ...% Fextx (n)
     1 1 1 1 ...% Fexty (n)
     1 1 1 1 ...% Mext (n)     
     0 0 0 0 ...% fidd (n)
     1 1 ];% basedd (2)
V=...
    [nan nan nan nan 0 ...% Frx (n+1)
     nan nan nan nan 0 ...% Fry (n+1)
     Mjo ...% M (n+1) % Mjo is [1x5]
     Fextx ...% Fextx (n)
     Fexty ...% Fexty (n)
     Mext  ...% Mext (n)
     nan nan nan nan ...% fidd (n)
     0 0 ];% basedd (2)

%constraint to keep xcomdd=0
% help variables:
% phi=phi(:)';
% phid=phid(:)';
% k(1)=-(m(1)*d(1)+m(2)*L(1)+m(3)*L(1));
% k(2)=-(m(2)*d(2)+m(3)*L(2));
% k(3)=-m(3)*d(3);
% aconstraints
% Acon=zeros(1,7*nseg+5);
% Acon(6*nseg+4:7*nseg+3)=k.*sin(phi);
% bconstraints
% bcon=sum(-k.*cos(phi).*phid.^2);
% update K and V; we make the hip-moment unknown to accomodate the
% constraint:
%K(3*nseg+2)=0;
%V(3*nseg+2)=nan;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% segdyn aanroepen
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[segdynstatedot, Vnew, succes] = segdyn(segdynstate, parms.segparms, K, V);
%[segdynstatedot, Vnew, succes] = segdyn(segdynstate, parms.segparms, K, V,Acon,bcon);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% applicatiespecifiek: statedot maken uit segdynstatedot
% in veel gevallen zullen segdynstatedot en statedot identiek zijn
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


statedot=[segdynstatedot(:)];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% applicatiespecifiek: uitgangen y berekenen
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% y is een KOLOMvector van willekeurige lengte die alle interessant geachte afhankelijke variabelen
% bevat (bijv reactiekrachten, versnellingen, energietermen, zwaartepuntpositie etc etc)

if parms.calculate_outputs
    output=[Vnew(:)];
end
