function [cm,cmp,cmdp,cmseg,cmdseg,cmddseg]=comtot(fitot,fiptot,fidptot,pos,posp,posdp,segparms)
%function [cm,cmp,cmdp]=comtot(fitot,fiptot,fidptot,pos,posp,posdp,segparms)
%
% INPUTS:	
%		fitot	containing angles, each sample should form one row
%       fiptot, fidptot time derivatives of fitot, may be empty
%		pos	position of base (provide emty argument if base is grounded at 0,0)
%       posp, posdp time derivatives of pos, may be empty
%       segparms  struct containing L,d,m,j
% OUTPUTS:	
%       cm,   nstep by 2 matrix containing x- and y-position of total cm
%       cmp,  nstep by 2 matrix containing x- and y-velocity of total cm
%       cmdp, nstep by 2 matrix containing x- and y-acceleration of total cm
% Note: number of outputs actually calculated depends on inputs provided
%			input fi                  => output: cm
%			input fi and fip          => output: cm and cmp
%			input fi and fip and fidp => output: cm and cmp and cmdp

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
    do_cm_velocity=1;
else
    do_cm_velocity=0;
end

if do_cm_velocity & all(size(fidptot)==[nstep nseg]) & all(size(posdp)==[nstep 2])
    do_cm_acceleration=1;
else
    do_cm_acceleration=0;
end

cm=zeros(nstep,2);
if do_cm_velocity
    cmp=zeros(nstep,2);
end

if do_cm_acceleration
    cmdp=zeros(nstep,2);
end

g=-9.81;
% doubling m makes calculations easier
mtot=sum(m);
m2=[m;m];

%?%? two rows, one column per segment
for i=1:nstep
    fi=fitot(i,:);fi=fi(:);
    dcfi=d.*cos(fi);dsfi=d.*sin(fi);
    lcfi=l.*cos(fi);lsfi=l.*sin(fi);
    
    % positie cm
    cmx(1,1)= dcfi(1);
    cmy(1,1)= dsfi(1);
    for j=2:nseg
        cmx(1,j)= sum(lcfi(1:j-1))+dcfi(j);
        cmy(1,j)= sum(lsfi(1:j-1))+dsfi(j); 
    end
    cm(i,1)=cmx*m/mtot;
    cm(i,2)=cmy*m/mtot;
    cmseg(i,:)=[cmx cmy];
    % velocities
    if do_cm_velocity
        fip=fiptot(i,:);fip=fip(:);
        % generate block representing first order transfer function 
        xblock = -[diag(dsfi)+tril(ones(nseg),-1)*diag(lsfi)];
        yblock =  [diag(dcfi)+tril(ones(nseg),-1)*diag(lcfi)];
        block  =  [xblock;yblock];
        % calculate cm velocity
        cmpseg=block*fip;
        cmp(i,1)=sum(cmpseg.*m2/mtot.*([ones(nseg,1);zeros(nseg,1)]));
        cmp(i,2)=sum(cmpseg.*m2/mtot.*([zeros(nseg,1);ones(nseg,1)]));
        cmdseg(i,:)=cmpseg'; % kkl
    end
    
    % accelerations
    if do_cm_acceleration
        fidp=fidptot(i,:);fidp=fidp(:);
        % generate block representing second order transfer function 
        xblock2 = [diag(dcfi)+tril(ones(nseg),-1)*diag(lcfi)];
        yblock2 = [diag(dsfi)+tril(ones(nseg),-1)*diag(lsfi)];
        block2  = [xblock2;yblock2];
        % calculate cm acceleration
        cmdpseg = block*fidp - block2*(fip.*fip);
        cmdp(i,1)=sum(cmdpseg.*m2/mtot.*([ones(nseg,1);zeros(nseg,1)]));
        cmdp(i,2)=sum(cmdpseg.*m2/mtot.*([zeros(nseg,1);ones(nseg,1)]));
        cmddseg(i,:)=cmdpseg';
    end
end

cm=cm+pos;
if do_cm_velocity, cmp=cmp+posp; end
if do_cm_acceleration, cmdp=cmdp+posdp; end

