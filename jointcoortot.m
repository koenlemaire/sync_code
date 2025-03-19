function [xpos,ypos, xposp, yposp,xposdp, yposdp]=jointcoortot(fitot,fiptot,fidptot,pos,posp,posdp,segparms)
%function [cm,cmp,cmdp]=comtot(fitot,fiptot,fidptot,pos,posp,posdp,segparms)
%
% INPUTS:	
%		fitot	containing angles, each sample should form one row
%       fiptot, fidptot time derivatives of fitot, may be empty
%		pos	position of base (provide emty argument if base is grounded at 0,0)
%       posp, posdp time derivatives of pos, may be empty
%       segparms  struct containing L,d,m,j
% OUTPUTS:	
%       xpos,    nstep by nseg+1 matrix containing x-coordinates of joints
%       ypos     nstep by nseg+1 matrix containing y-coordinates of joints
%       xposp    nstep by nseg+1 matrix containing x-velocities of joints
%       yposp    nstep by nseg+1 matrix containing y-velocities of joints
%       xposdp   nstep by nseg+1 matrix containing x-accelerations of joints
%       yposdp   nstep by nseg+1 matrix containing y-accelerations of joints
%
% Note: number of outputs actually calculated depends on inputs provided

[nstep,nseg]=size(fitot);

l=segparms.L(:);
d=segparms.d(:);
m=segparms.m(:);
j=segparms.j(:);

if length(l)~=nseg   
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

if isempty(pos)
    disp('assuming grounded base');
    pos=zeros(nstep,2);posp=zeros(nstep,2);posdp=zeros(nstep,2);
end

if all(size(fiptot)==[nstep nseg]) & all(size(posp)==[nstep 2])
    do_velocity=1;
else
    do_velocity=0;
end

if do_velocity & all(size(fidptot)==[nstep nseg]) & all(size(posdp)==[nstep 2])
    do_acceleration=1;
else
    do_acceleration=0;
end

xpos=zeros(nstep,nseg+1);
ypos=zeros(nstep,nseg+1);
if do_velocity
    xposp=zeros(nstep,nseg+1);
    yposp=zeros(nstep,nseg+1);
end

if do_acceleration
    xposdp=zeros(nstep,nseg+1);
    yposdp=zeros(nstep,nseg+1);
end

%?%? two rows, one column per segment
for i=1:nstep
    fi=fitot(i,:);fi=fi(:);
    lcfi=l.*cos(fi);lsfi=l.*sin(fi);
    
    % positie 
    xpos(i,1)= pos(i,1);
    ypos(i,1)= pos(i,2);
    for j=1:nseg
        xpos(i,j+1)= xpos(i,j)+lcfi(j);
        ypos(i,j+1)= ypos(i,j)+lsfi(j);
    end
    
    % velocities
    if do_velocity
        xposp(i,1)= posp(i,1);
        yposp(i,1)= posp(i,2);
        for j=1:nseg
            xposp(i,j+1)= xposp(i,j)-lsfi(j)*fiptot(i,j);
            yposp(i,j+1)= yposp(i,j)+lcfi(j)*fiptot(i,j);
        end
    end
    
    % accelerations
    if do_acceleration
        xposdp(i,1)= posdp(i,1);
        yposdp(i,1)= posdp(i,2);
        for j=1:nseg
            xposdp(i,j+1)= xposdp(i,j)-lsfi(j)*fidptot(i,j)-lcfi(j)*fiptot(i,j)^2;
            yposdp(i,j+1)= yposdp(i,j)+lcfi(j)*fidptot(i,j)-lsfi(j)*fiptot(i,j)^2;
        end
    end
end
