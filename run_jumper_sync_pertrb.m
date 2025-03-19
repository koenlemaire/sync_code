%% forward -- parameters
% the model is an inverted double pendulum, with segment masses and
% inertias close to human legs and trunk. The pendulum is kept upright by
% torsional spring-dampers at the ankle and hip/knee, with parameters such
% that the system is subcritically damped, and barely stable. The
% initial condition is .1 rad out of the equilibrium position for both
% segments.

clear all;
close all; clc 

% control switches:
animate=0;
Ncolor=[0 0 0];
alphaMod=linspace(0.4,1,3);alphaMod=fliplr(alphaMod);
lblx=0.05;
lbly=1.05;

% based on simulations made by Richard Casius / Knoek van Soest
tmp1=load('segpar');
tmp2=load('depvarj');
tmp3=load('state');

% segmentparameters
parms.segparms.L=tmp1.len;
parms.segparms.d=tmp1.d;
parms.segparms.m=tmp1.mass;
parms.segparms.j=tmp1.I;

L=parms.segparms.L;
m=parms.segparms.m;
d=parms.segparms.d;

nseg=length(parms.segparms.L);

% COM height at straight up standing 1.0109 m
standing_com=1.0409; % [m]

% normalization constants
force_norm=sum(parms.segparms.m)*9.81; % [N] force normalization constant
length_norm=sum(parms.segparms.L(2:3)); % [m] length normalization constant, lower leg + upper leg length
time_norm=sqrt(parms.segparms.L(1)/9.81); % [s] time normalization constant

moment_norm=force_norm*length_norm; % [Nm]
power_norm=force_norm*length_norm/time_norm; % [W]
work_norm=moment_norm;

%% controls
parms.controls=tmp2.joint_depvar(11:15,1:end-20); % joint moments
parms.t_control=tmp3.time(1:end-20);

% add 0.15 s period of constant control in the beginning and end : FIRST
% CALCULATE STATIC CONDITION
%parms.controls=[repmat(parms.controls(:,1),1,15) parms.controls];
%parms.t_control=[linspace(0,.14,15) parms.t_control+0.15];


% initial conditions
phi0=tmp3.state(1:4,1);
phid0=tmp3.state(5:8,1);
base0=tmp3.state(21:22,1);
based0=tmp3.state(23:24,1);

state0=[phi0; phid0; base0; based0];

%% forward -- take-off phase
% first calculate equilibrium condition of quiet standing:
phidd0=phid0; % zeros
basedd0=based0; % zeros
[y0]=inv_jumper_equilibrium(0,state0,phidd0,basedd0,parms);
Mjo_eq=y0(2*nseg+3:3*nseg+3); % equilibrium joint moments

% update control to start with equilibrium joint moments:
parms.controls=[repmat(Mjo_eq,1,16) parms.controls];
%parms.t_control=[linspace(0,.05,51) parms.t_control+0.0501];
parms.t_control=[linspace(0,.015,16) parms.t_control+0.0151];

parms.calculate_outputs=true;
[stated0,y0]=frwd_jumper(0,state0,parms);

tspan=[0 5];
odeopt=odeset('abstol',1e-10,'reltol',1e-10,'events',@jumper_event);

tic
[t_to,state_to]=ode113(@frwd_jumper,tspan,state0,odeopt,parms);
toc
for i=1:length(t_to)
    [stated_to(:,i),y_to(:,i)]=frwd_jumper(t_to(i),state_to(i,:),parms);
end
take_off=t_to(end); % [s] take-off time 
stated_to=stated_to';y_to=y_to';
%% forward -- flight phase
tspan=[t_to(end) t_to(end)+.015];
odeopt=odeset('abstol',1e-10,'reltol',1e-10,'events',@jumper_event_flight);

tic
[t_flight,state_flight]=ode113(@frwd_jumper_flight,tspan,state_to(end,:),odeopt,parms);
toc
for i=1:length(t_flight)
    [stated_flight(:,i),y_flight(:,i)]=frwd_jumper_flight(t_to(i),state_flight(i,:),parms);
end
stated_flight=stated_flight';y_flight=y_flight';
%% forward -- dependent variables
% merging pre and post flight: 
t = [t_to; t_flight(2:end)];
state = [state_to; state_flight(2:end,:)];
stated = [stated_to; stated_flight(2:end,:)];
y = [y_to; y_flight(2:end,:)];

% interpolate to equidistant time base:
dt=.0005; % 1000 Hz sampling rate
t_new=0:dt:t(end);
for i=1:length(state0)
    state_new(:,i)=interp1(t,state(:,i),t_new);
    stated_new(:,i)=interp1(t,stated(:,i),t_new);
end
for i=1:length(y0)
    y_new(:,i)=interp1(t,y(:,i),t_new);
end
t=t_new;
state=state_new;
stated=stated_new;
y=y_new;



% dependent variables:
grf=[y(:,1) y(:,nseg+2)]; % ground reaction force
Fjo=[y(:,2:nseg) y(:,nseg+3:2*nseg+1)]; % joint reaction forces
Fres=-[y(:,nseg+1) y(:,2*nseg+2)]; % residual force (end of chain, acting on )
Mjo=y(:,2*nseg+3:3*nseg+3); % joint moments
Mres=-y(:,3*nseg+3); % residual moment (end of chain)
grm=y(:,2*nseg+3); % moment from ground on subject

% unraveling state:
phi=state(:,1:nseg);
phid=state(:,nseg+1:2*nseg);
phidd=stated(:,nseg+1:2*nseg);
base=state(:,2*nseg+1:2*nseg+2);
based=state(:,2*nseg+3:2*nseg+4);
basedd=stated(:,2*nseg+3:2*nseg+4);

% joint angles, angular velocity
phiJo=phi(:,2:nseg)-phi(:,1:nseg-1);
phidJo=phid(:,2:nseg)-phid(:,1:nseg-1);

% total and segment com kinematics (based on kinematics ...)
[com,comd,comdd,cmseg,cmdseg,cmddseg]=comtot(phi,phid,phidd,base,based,basedd,parms.segparms);
% cmseg = [x x y y]
% joint kinematics
[xpos,ypos,xposd,yposd,xposdd,yposdd]=jointcoortot(phi,phid,phidd,base,based,basedd,parms.segparms);

% energy terms
[epot,ekinx,ekiny,erot,etot]=energy(phi,phid,base,based,parms.segparms);
etot=etot-etot(1);

% joint work [not including residual force ...]
for i=1:nseg-1
    Pjo(:,i)=Mjo(:,i+1).*(phid(:,i+1)-phid(:,i));
    Wjo(:,i)=cumtrapz(phi(:,i+1)-phi(:,i),Mjo(:,i+1));
end
Pjo_sum=sum(Pjo,2);
Wjo_sum=sum(Wjo,2);

% residual work
W_Fres=cumtrapz(xpos(:,end),Fres(:,1)) + cumtrapz(ypos(:,end),Fres(:,2));
W_Mres=cumtrapz(phi(:,end),Mres);
Wres=W_Fres+W_Mres;

% COM work
acom=grf/sum(parms.segparms.m);
acom(:,2)=acom(:,2)-9.81;
vcom=[cumtrapz(t,acom(:,1)) cumtrapz(t,acom(:,2))];
Pcom=sum(grf.*vcom,2); % [W] COM power, using kinetics based COM
Wcom=cumtrapz(t,Pcom); % [J] COM work

% Peripheral work DO MINUS VCOM RATHER THAN COMD???
for i=1:nseg
    Wper_transx(:,i)=.5*parms.segparms.m(i)*(cmdseg(:,i)-vcom(:,1)).^2;
    Wper_transy(:,i)=.5*parms.segparms.m(i)*(cmdseg(:,i+nseg)-vcom(:,2)).^2;
end
Wper_tot=sum(Wper_transx,2)+sum(Wper_transy,2)+sum(erot,2);

% COM + peripheral work
WcomPLUSper=Wcom+Wper_tot;

% Art soft tissue work
Wsoft=WcomPLUSper-Wjo_sum;

% real soft tissue work
Wsoft_real=cumtrapz(base(:,1),grf(:,1))+cumtrapz(base(:,2),grf(:,2));

% residual work check
Wres_check=etot-Wjo_sum-Wsoft_real;

% calculate jump height:
itmp=find(comdd(:,2)<-9.8);
itmp=itmp(end);
com_to=com(itmp,2);
vcom_to=comd(itmp,2);
com_apex=com_to+(.5*vcom_to^2)/9.81;
jump_height=com_apex-standing_com

%% forward figures
figure(1);subplot(211);plot(t-take_off,phi)
title('angles')
xlabel('Time [s]')
ylabel('Angle [rad]')

% phi, phid
subplot(212);plot(phi,phid)
title('phase diagrams of segments')
xlabel('Angle [rad]')
ylabel('Angular velocity [rad/s]')

figure
subtmp=subplot(3,2,[1 3 5]);
yoff=.4;
%tmp.Position(2)=tmp.Position(2)+.1;
%tmp.Position(4)=tmp.Position(4)-.1;
%plot(xpos(1,:),ypos(1,:),'k-'); hold on
plot(xpos(450,:)-xpos(1,1),ypos(450,:)+yoff,'k-','linewidth',1); hold on % stick figure
circr=.1;
xcirc=xpos(450,end)-xpos(1,1)+cos(phi(450,end))*circr;
ycirc=ypos(450,end)+yoff+sin(phi(450,end))*circr;
rectangle('Curvature',[1 1],'Position',[xcirc-circr ycirc-circr 2*circr 2*circr],'linewidth',1);

plot([-.3 .3],[yoff yoff],'color',[0 0 0 .8],'linewidth',.7);hold on
nline=12;
tmp=linspace(-.3,.3,nline);
for iLine=1:nline
    plot([tmp(iLine)-.05 tmp(iLine)],[yoff-.05 yoff],'color',[0 0 0 .8],'linewidth',.7)
end
plot([0 0],[0 yoff-.01],'r','linewidth',1)
plot([0 -.025],[yoff-.01 yoff-.06],'r','linewidth',1)
plot([0 +.025],[yoff-.01 yoff-.06],'r','linewidth',1)
text(subtmp,0.2,1.02,'A','units','normalized')
text(subtmp,1.35,1.02,'B','units','normalized')
text(subtmp,1.35,.66,'C','units','normalized')

%plot(,'r','linewidth',1)
set(gca, 'Xcolor', 'w', 'Ycolor', 'w')
set(gca, 'XTick', []);
set(gca, 'YTick', []);
box off        
%plot(xpos(590,:),ypos(590,:),'k-'); hold on
axis equal


subplot(3,2,2)
idxtmp=t>.015 & t<t(end)-.015;
plot(t(idxtmp)-take_off,phiJo(idxtmp,1),'k-')
xlim([t(1) t(end)]-take_off)
xlabel('Time (s)')
ylabel('Ankle angle (rad)')
box off
%text(lblx,lbly,'B','units','normalized')

tmp=subplot(3,2,[4,6]);
%tmp=axes;
yyaxis('left')
tmp.YColor=[0 0 0];
plot(t-take_off,grf(:,2)/force_norm,'-','color',[Ncolor alphaMod(2)]); hold on
plot(t-take_off-.015,grf(:,2)/force_norm,'-','color',[Ncolor alphaMod(3)]); hold on
plot(t-take_off+0.015,grf(:,2)/force_norm,'-','color',[Ncolor alphaMod(1)]); hold on
text(-.16,2.6,'$+15 \, \mathrm{ms}$','fontsize',8)
text(-.29,2.6,'$-15 \, \mathrm{ms}$','fontsize',8)

xlim([t(1) t(end)]-take_off)
ylim([-0.5 3])
xlabel('Time (s)')
ylabel('Force ($Mg$)')
%ltmp=legend('horizontal','vertical');
%ltmp.AutoUpdate='off';
%plot(t(1)-take_off+[0.015 0.015],[-.5 3],'k--',[0 0],[-.5 3],'k--')
ptmp=patch([t(1)-take_off+[0.015 0.015] 0 0],[-.49 3 3 -.49],[.7 .7 .7]);
ptmp.FaceAlpha=.2;
ptmp.EdgeColor=[1 1 1];
ptmp.EdgeAlpha=[0];
yyaxis right
tmp.YColor=[0 0 0];
tmp.YLim=[-0.5 3]*force_norm;
ylabel('(N)')
box off


if animate
    figure(99)
    pause
    for i =1:10:length(t)
        % ceiling
        plot([-1 1],[0 0],'k','linewidth',1);hold on
        nline=12;
        tmp=linspace(-1,1,nline);
        for iLine=1:nline
            plot([tmp(iLine) tmp(iLine)-.1],[0 -.1],'k','linewidth',1)
        end
        % kinematics
        plot(xpos(i,:),ypos(i,:),'k','linewidth',1.5); hold off
        %plot(inv_xpos(i,:),inv_ypos(i,:),'r--','linewidth',1.5); hold off
        ylim([-.5,2])
        axis equal
        drawnow
    end
end

%% time lag perturbation to kinematic variables
% The kinematics are made to lag the force recording, because in this
% scenario the timing of the analysis window will be determined by the
% forces, ie the force is the leading signal
lag=.015; % [s] time lag range
delayVec=linspace(-lag,lag,31);
for iLag=1:length(delayVec)
    % bottom steps to ensure wobbly mass is correctly implemented ...
    inv_parms=parms;
    % apply any perturbations to the measurement (to do apply perturbation to the
    % model, change the forward part of this program ....)
    delay=-delayVec(iLag); % nr of samples force is lagging mocap [ms]
    % perturbations to parameters
    %inv_parms.segparms.m(1)=inv_parms.segparms.m(1)-2;
    %inv_parms.segparms.m(2)=inv_parms.segparms.m(2)+2;
    idx_old=t>.0151 & t<t_to(end);
    t_old=t(idx_old);
    t_new=t_old+delay; % new time base
    
    plot_idx=max([1,delay]):min([length(t_new),length(t_new)+delay]);
    
    % cut away excess pieces
    %t_new(t_new<0)=[];
    %t_new(t_new>t(end))=[];
        
    % time lag perturbation to kinematics
    %inv_t=t(inv_idx);
    for i=1:length(state0)
        inv_state(:,i)=interp1(t,state(:,i),t_new);
        inv_stated(:,i)=interp1(t,stated(:,i),t_new);
    end
    
    % unravel new inverse state
    inv_phi=inv_state(:,1:nseg);
    inv_phid=inv_state(:,nseg+1:2*nseg);
    inv_phidd=inv_stated(:,nseg+1:2*nseg);
    inv_base=inv_state(:,2*nseg+1:2*nseg+2);
    inv_based=inv_state(:,2*nseg+3:2*nseg+4);
    inv_basedd=inv_stated(:,2*nseg+3:2*nseg+4);    
    
    % joint angles, angular velocity
    inv_phiJo=inv_phi(:,2:nseg)-inv_phi(:,1:nseg-1);
    inv_phidJo=inv_phid(:,2:nseg)-inv_phid(:,1:nseg-1);
    
    % perturbations to ground contact:
    inv_grf=grf(idx_old,:);
    inv_grm=grm(idx_old,:);
    
    % visualisation of perturbations:
    delay_col=[0 0 0 .7];
    if iLag==1 || iLag==length(delayVec)
        figure(3)
        subplot(211)
        plot(t-take_off,phi(:,2)-phi(:,1)+pi,'k');hold on
        plot(t_old-take_off,inv_phi(:,2)-inv_phi(:,1)+pi,'color',delay_col) % ground truth
        title('ankle angle')
        xlabel('Time [s]')
        ylabel('rad')
        legend('forward','inverse')
        
        subplot(212)
        plot(t-take_off,phid(:,2)-phid(:,1),'k');hold on
        plot(t_old-take_off,inv_phid(:,2)-inv_phid(:,1),'color',delay_col) % ground truth
        title('ankle angular velocity')
        xlabel('Time [s]')
        ylabel('rad/s')
        legend('forward','inverse')
    end
    %% inverse analysis    
    for i=1:length(t_new)
        [inv_y(:,i)]=inv_jumper(t_new(i),inv_state(i,:),inv_phidd(i,:),inv_basedd(i,:),inv_grf(i,:),inv_grm(i),inv_parms);
    end
    inv_y=inv_y';
    inv_Fjo=[inv_y(:,2:nseg) inv_y(:,nseg+3:2*nseg+1)]; % joint reaction forces
    inv_Fres=-[inv_y(:,nseg+1) inv_y(:,2*nseg+2)]; % residual force (end of chain)
    inv_Mjo=inv_y(:,2*nseg+3:3*nseg+3); % joint moments
    inv_y=inv_y';
    % total and segment com kinematics (based on kinematics ...)
    [inv_com,inv_comd,inv_comdd,inv_cmseg,inv_cmdseg,inv_cmddseg]=comtot(inv_phi,inv_phid,inv_phidd,inv_base,inv_based,inv_basedd,inv_parms.segparms);
    
    % joint kinematics
    [inv_xpos,inv_ypos,inv_xposd,inv_yposd,inv_xposdd,inv_yposdd]=jointcoortot(inv_phi,inv_phid,inv_phidd,inv_base,inv_based,inv_basedd,inv_parms.segparms);
    
    % residual moment about COM
    r_end_wrt_com=[inv_xpos(:,end) inv_ypos(:,end)]-inv_com;
    M_Fres_wrt_com=r_end_wrt_com(:,1).*inv_Fres(:,2) - r_end_wrt_com(:,2).*inv_Fres(:,1);
    inv_Mres=-inv_Mjo(:,end)+M_Fres_wrt_com; % residual moment about COM!
        
    % energy terms
    [inv_epot,inv_ekinx,inv_ekiny,inv_erot,inv_etot]=energy(inv_phi,inv_phid,inv_base,inv_based,inv_parms.segparms);
    inv_etot=inv_etot-inv_etot(1);
    
    % joint work [not including residual force ...]
    for i=1:nseg-1
        inv_Pjo(:,i)=inv_Mjo(:,i+1).*(inv_phid(:,i+1)-inv_phid(:,i));
        inv_Wjo(:,i)=cumtrapz(inv_phi(:,i+1)-inv_phi(:,i),inv_Mjo(:,i+1));
    end
    inv_Wjo_net(iLag,:)=inv_Wjo(end,:); % net joint work
    inv_Pjo_sum=sum(inv_Pjo,2);
    inv_Wjo_sum=sum(inv_Wjo,2);
    
    % residual work
    inv_W_Fres=cumtrapz(inv_xpos(:,end),inv_Fres(:,1)) + cumtrapz(inv_ypos(:,end),inv_Fres(:,2));
    inv_W_Mres=cumtrapz(inv_phi(:,end),inv_Mres);
    inv_Wres=inv_W_Fres+inv_W_Mres;
    inv_Wres_check=inv_etot-inv_Wjo_sum;
    
    % COM work based on grf (assuming vcom(0)=0)
    inv_acom=inv_grf/sum(inv_parms.segparms.m); % [m/s^2]
    inv_acom(:,2)=inv_acom(:,2)-9.81; % [m/s^2]
    inv_vcom=[cumtrapz(t_old,inv_acom(:,1)) cumtrapz(t_old,inv_acom(:,2))]; % [m/s]
    %inv_vcom=[cumtrapz(t,inv_acom(:,1))+vcom(n_sample+1,1) cumtrapz(t,inv_acom(:,2))+vcom(n_sample+1,2)]; % [m/s]
    %inv_vcom=[cumtrapz(t,inv_acom(:,1))+comd(n_sample+1,1) cumtrapz(t,inv_acom(:,2))+comd(n_sample+1,2)]; % [m/s]
    inv_Pcom=sum(inv_grf.*inv_vcom,2); % [W] COM power
    inv_Wcom=cumtrapz(t_old,inv_Pcom); % [J] COM work
    
    % Peripheral work
    for i=1:nseg
        inv_Wper_transx(:,i)=.5*inv_parms.segparms.m(i)*(inv_cmdseg(:,i)-inv_vcom(:,1)).^2;
        inv_Wper_transy(:,i)=.5*inv_parms.segparms.m(i)*(inv_cmdseg(:,i+nseg)-inv_vcom(:,2)).^2;
    end
    inv_Wper_tot=sum(inv_Wper_transx,2)+sum(inv_Wper_transy,2)+sum(inv_erot,2);
    
    % COM + peripheral work
    inv_WcomPLUSper=inv_Wcom+inv_Wper_tot;
    
    % soft tissue work
    inv_Wsoft=inv_WcomPLUSper-inv_Wjo_sum;   
    
    %% joint angle, torque and power error
    % surprisingly large change in joint torques; likely depends on
    % frequency of the signal, but that is not too outrageous in this
    % example ...
    % put joint work error in
    angle_error_tot(iLag,:)=sqrt(mean((phiJo(idx_old,:)-inv_phiJo(:,:)).^2));
    
    Mjo_error(iLag)=sqrt(mean((Mjo(idx_old,2)-inv_Mjo(:,2)).^2));
    Mjo_error_tot(iLag,:)=sqrt(mean((Mjo(idx_old,:)-inv_Mjo(:,:)).^2))/moment_norm;

    % similar to joint torque, but now the errors accumulate over time; we
    % have taken a typical (step) time to analyse here..
    Pjo_error(iLag)=sqrt(mean((Pjo(idx_old,1)-inv_Pjo(:,1)).^2));
    Pjo_error_tot(iLag,:)=sqrt(mean((Pjo(idx_old,:)-inv_Pjo(:,:)).^2))/power_norm;
    
    Fres_mag(iLag)=mean(sqrt(sum(inv_Fres.^2,2))); % [N] mean of Fres magnitude
    Mres_mag(iLag)=mean(abs(inv_Mres));
    
    %% 3x3
    if  iLag==1 || iLag==length(delayVec) || delayVec(iLag) == 0 
        figure(49)
        joNames=[{'Ankle'},{'Knee'},{'Hip'}];
        for iJo=1:3
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % angle
            tmp=subplot(3,3,iJo);
            yyaxis left
            tmp.YColor=[0 0 0];
            if iLag == 1                
                %plot(t_old-take_off,inv_phiJo(:,iJo),'-','color',[Ncolor alphaMod(1)]); hold on % ground truth
            elseif iLag==length(delayVec)
                %plot(t_old-take_off,inv_phiJo(:,iJo),'-','color',[Ncolor alphaMod(3)]) % ground truth                
                if iJo==1
                    %annotation('textarrow',[0.20 0.14],[0.87 0.87],'String','$+ \Delta t \,$')
                end
            else         
                plot(t-take_off,phiJo(:,iJo),'-','color',[Ncolor]);hold on
            end
            title(joNames(iJo))
            if iJo==1
                ylabel('Joint angle (rad)')
            end
            tmp2=tmp.YLim;
            yyaxis right
            tmp.YColor=[0 0 0];
            tmp.YLim=tmp2*180/pi;
            box off
            if iJo==3
                ylabel('(deg)')
            end
            

            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % moment
            tmp=subplot(3,3,iJo+3);
            yyaxis left
            tmp.YColor=[0 0 0];
            if iLag == 1                            
                plot(t_old-take_off,inv_Mjo(:,iJo+1)/moment_norm,'-','color',[Ncolor alphaMod(1)]); hold on % inverse result                
            elseif iLag == length(delayVec)
                plot(t_old-take_off,inv_Mjo(:,iJo+1)/moment_norm,'-','color',[Ncolor alphaMod(3)]) % inverse result                
            else
                plot(t-take_off,Mjo(:,iJo+1)/moment_norm,'-','color',[Ncolor alphaMod(2)]); hold on % ground truth
            end
            box off
            %tmp=text(1,1,'$\Delta t=15 \mathrm{ms}$','interpreter','latex','units','centimeters');
            %tmp=text(1,2,'$\Delta t=-15 \mathrm{ms}$','interpreter','latex','units','centimeters');
            if iJo==1
                ylabel('Moment ($MgL$)')
            end
            xlim([t_old(1)-take_off 0])
            tmp2=tmp.YLim;
            yyaxis right
            tmp.YColor=[0 0 0];
            tmp.YLim=tmp2*moment_norm;
            box off
            if iJo==3
                ylabel('(Nm)')
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % power
            tmp=subplot(3,3,iJo+6);
            yyaxis left
            tmp.YColor=[0 0 0];
            if iLag == 1                            
                plot(t_old-take_off,inv_Pjo(:,iJo)/power_norm,'-','color',[Ncolor alphaMod(1)]); hold on % inverse result
            elseif iLag == length(delayVec)
                plot(t_old-take_off,inv_Pjo(:,iJo)/power_norm,'-','color',[Ncolor alphaMod(3)]) % inverse result
            else
                plot(t-take_off,Pjo(:,iJo)/power_norm,'-','color',[Ncolor alphaMod(2)]); hold on % ground truth
            end
            %tmp=text(1,1,'$\Delta t=15 \mathrm{ms}$','interpreter','latex','units','centimeters');
            %tmp=text(1,2,'$\Delta t=-15 \mathrm{ms}$','interpreter','latex','units','centimeters');
            box off
            if iJo==1
                ylabel('Power ($Mg \sqrt{gL}$)')
            end            
            xlabel('Time (s)')
            xlim([t_old(1)-take_off 0])
            tmp2=tmp.YLim;
            yyaxis right
            tmp.YColor=[0 0 0];
            tmp.YLim=tmp2*power_norm;
            box off
            if iJo==3
                ylabel('(W)')
            end
            
        end
        %ylabel('Moment [m*g*LL]')
        %legend('forward','inverse')
    end
    
    
end
%% average errors
% 3x3
figure(50)

tmp=subplot(221);
yyaxis left
tmp.YColor=[0 0 0];
plot(delayVec*1000,angle_error_tot,'k')
%xlabel('delay [ms]')
ylabel('Joint angle (rad)')
text(tmp,lblx,lbly,'A','units','normalized')
%title('A')
legend('ankle','knee','hip','Location','north')
legend box off
box off
xlim([-15 15])
xticks([-15:5:15])
tmp2=tmp.YLim;
yyaxis right
ylabel('(deg)')
tmp.YColor=[0 0 0];
tmp.YLim=tmp2*180/pi;
box off

tmp=subplot(222);
yyaxis left
tmp.YColor=[0 0 0];
plot(delayVec*1000,Mjo_error_tot(:,2:4),'k')
%xlabel('delay [ms]')
ylabel('Moment ($MgL$)')
text(tmp,lblx,lbly,'B','units','normalized')
%legend('ankle','knee','hip')
xlim([-15 15])
xticks([-15:5:15])
tmp2=tmp.YLim;
yyaxis right
ylabel('(Nm)')
tmp.YColor=[0 0 0];
tmp.YLim=tmp2*moment_norm;
box off

tmp=subplot(223);
yyaxis left
tmp.YColor=[0 0 0];
plot(delayVec*1000,Pjo_error_tot,'k')
xlabel('Delay (ms)')
ylabel('Power ($Mg \sqrt{gL}$)')
text(tmp,lblx,lbly,'C','units','normalized')
xlim([-15 15])
xticks([-15:5:15])
%legend('ankle','knee','hip')
tmp2=tmp.YLim;
yyaxis right
ylabel('(W)')
tmp.YColor=[0 0 0];
tmp.YLim=tmp2*power_norm;
box off


%print('3x3','-dsvg')
%print('3x3','-dtiff')

% work errors

tmp=subplot(224);
yyaxis('left')
tmp.YColor=[0 0 0];
plot(delayVec*1000,inv_Wjo_net/moment_norm,'k') % ground truth
%title('Net joint work')
xlabel('Delay (ms)')
%ylabel('Net work [$m \cdot g \cdot \mathrm{LL}$]')
ylabel('Net work ($MgL$)')
text(tmp,lblx,lbly,'D','units','normalized')
xlim([-15 15])
xticks([-15:5:15])
%legend('ankle','knee','hip','Location','best')
ymax=.4;
ylim([0 ymax])
yyaxis right
ylabel('(J)')
ylim([0 ymax]*moment_norm)
tmp.YColor=[0 0 0];
box off
%xlim([-lag lag])

figure(70)
tmp=subplot(121);
yyaxis('left')
tmp.YColor=[0 0 0];
plot(delayVec*1000,Fres_mag/force_norm,'k') % ground truth
%title('Net joint work')
xlabel('Delay (ms)')
%ylabel('Net work [$m \cdot g \cdot \mathrm{LL}$]')
ylabel('Residual force ($Mg$)')
text(tmp,lblx,lbly,'A','units','normalized')
xlim([-15 15])
xticks([-15:5:15])
%legend('ankle','knee','hip','Location','best')
ymax=.25;
ylim([0 ymax])
yyaxis right
ylabel('(N)')
ylim([0 ymax]*force_norm)
tmp.YColor=[0 0 0];
box off

tmp=subplot(122);
yyaxis('left')
tmp.YColor=[0 0 0];
plot(delayVec*1000,Mres_mag/moment_norm,'k') % ground truth
%title('Net joint work')
xlabel('Delay (ms)')
%ylabel('Net work [$m \cdot g \cdot \mathrm{LL}$]')
ylabel('Residual moment ($MgL$)')
text(tmp,lblx,lbly,'B','units','normalized')
xlim([-15 15])
xticks([-15:5:15])
%legend('ankle','knee','hip','Location','best')
ymax=.05;
ylim([0 ymax])
yyaxis right
ylabel('(N m)')
ylim([0 ymax]*moment_norm)
tmp.YColor=[0 0 0];
box off

% plan: make figures of real vs inverse things for all interesting
% variables, do this for different types of perturbations.

% do:
% i) drop bag of rice
% ii) any rigid motion
% iii) soft tissue in attachment and/or soft tissue in
% iv) make kinematic perturbation physical; skin marker motion via spring damper
% misestimates of joint centre (low prio)
% v) test effect of delay of force data wrt phasespace data ...

% look at:
% soft tissue work
% joint torques
% residual force