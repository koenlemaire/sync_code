function [epot,ekinx,ekiny,erot,etot]=energy(fitot,fiptot,pos,posp,segparms)
%function [epot,ekinx,ekiny,erot,etot]=energy(fitot,fiptot,pos,posp,parms)
%
% INPUTS
% 
% fitot     matrix of segment angles, each sample should form one row
% fiptot    matrix of segment angular velocities, each sample should form one row
% pos       matrix of base position (x y), each sample should form one row 
%           (base is assumed to be fixed at (0,0) if pos is not supplied)
% posp      matrix of base velocity (x y), each sample should form one row
%           (base velocity is assumed to be (0,0) if posp is not supplied)
% segparms  STRUCT containing segment parameter vectors L,d,m,j
%
% OUTPUTS
% epot      matrix of segment potential energies (assuming positive y is antigravity)
% ekinx     matrix of segment com horizontal kinetic energies
% ekiny     matrix of segment com vertical kinetic energies
% erot      matrix of segment rotational energies
% etot      vector of total mechanical energy

[nr,nc]=size(fitot);nstep=nr;nseg=nc;

if any(size(fiptot)~=size(fitot))
    disp('ERROR: dimensions of fi and fip should be equal')
    return
end

if isempty(pos)|isempty(posp)
    disp('assuming grounded base');
    pos=zeros(nstep,2);posp=zeros(nstep,2);
else
    if any(size(pos)~=[nstep 2])    
        disp('ERROR: dimensions of input arg 3 are inconsistent')
        return
    end
    
    if any(size(posp)~=[nstep 2])    
        disp('ERROR: dimensions of input arg 4 are inconsistent')
        return
    end
end

L=segparms.L(:);
d=segparms.d(:);
m=segparms.m(:);
j=segparms.j(:);

if length(L)~=nseg   
    disp('ERROR: dimension of parms.L is incorrect')
    return
end

if length(d)~=nseg   
    disp('ERROR: dimension of parms.d is incorrect')
    return
end
if length(m)~=nseg   
    disp('ERROR: dimension of parms.m is incorrect')
    return
end
if length(j)~=nseg   
    disp('ERROR: dimension of parms.j is incorrect')
    return
end

g=9.81;
% doubling m makes calculations easier
mtot=sum(m);
m2=[m;m];

% initializing output matrices
epot=zeros(nstep,nseg);
ekinx=zeros(nstep,nseg);
ekiny=zeros(nstep,nseg);
erot=zeros(nstep,nseg);
etot=zeros(nstep,1);

%?%? two rows, one column per segment
for i=1:nstep
    fi=fitot(i,:);fi=fi(:);
    dcfi=d.*cos(fi);dsfi=d.*sin(fi);
    lcfi=L.*cos(fi);lsfi=L.*sin(fi);
    
    % positie cm
    epot(i,1)= (pos(i,2)+dsfi(1))*m(1)*g;
    for ii=2:nseg
        epot(i,ii)= (pos(i,2)+sum(lsfi(1:ii-1))+dsfi(ii))*m(ii)*g; 
    end
    
    % velocities
    fip=fiptot(i,:);fip=fip(:);
    % generate block representing first order transfer function 
    xblock = -[diag(dsfi)+tril(ones(nseg),-1)*diag(lsfi)];
    yblock =  [diag(dcfi)+tril(ones(nseg),-1)*diag(lcfi)];
    block  =  [xblock;yblock];
    % calculate cm velocity
    cmpseg=block*fip;
    ekinx(i,:)=(0.5*m.*(posp(i,1)+cmpseg(1:nseg)).^2)';
    ekiny(i,:)=(0.5*m.*(posp(i,2)+cmpseg(nseg+1:2*nseg)).^2)';
    erot(i,:)=(0.5*j.*fip.^2)';
end

if nseg==1
    etot=epot+ekinx+ekiny+erot;
else
    etot=sum(epot'+ekinx'+ekiny'+erot')';
end
