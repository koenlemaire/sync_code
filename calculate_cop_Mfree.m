function [copl, copr, Mfree_l, Mfree_r] = calculate_cop_Mfree (grfl,grfr,grml,grmr,LFT_rc,RFT_rc,doFig)
% calculates COP and free moment based on grm and grf data and implements a
% correction based on foot trajectory in space. Assumes GRM and GRF are
% already transported to the MOCAP coordinate system!! 
copl=[-grml(:,2)./grfl(:,3) grml(:,1)./grfl(:,3) zeros(length(grfl),1)];
copr=[-grmr(:,2)./grfr(:,3) grmr(:,1)./grfr(:,3) zeros(length(grfr),1)];

if doFig
    fig=figure;
    fig.Name=['cop trajectory check'];
    subplot(221);plot(copl(:,1),copl(:,2),'.'); hold on
    plot(copr(:,1),copr(:,2),'.');
    title('uncorrected COP')
end
% COP correction:
% we are going to constrain the cop trajectory such that the cop is
% never more than foot_radius away from the foot COM. If the cop
% location is outside of foot_radius it is projected onto the circle
% edge. This ensures that the cop trajectory is continuous.

foot_radius=0.15; % [m] tolerable distance from com foot
%copl_raw=copl;
%copr_raw=copr;

% left foot
cop_vec=copl(:,1:2)-LFT_rc(:,1:2); % vector from com foot to cop
cop_dist=sqrt(sum(cop_vec.^2,2)); % distance of cop from com foot
cop_evec=[cop_vec(:,1)./cop_dist cop_vec(:,2)./cop_dist];
idx=cop_dist>foot_radius; % [m] Pythagoras]
copl(idx,1:2)=LFT_rc(idx,1:2) + foot_radius*cop_evec(idx,:);

% symmetry correction???? KKL: this is weird currently
copl(:,2)=copl(:,2)+0.00;
copr(:,2)=copr(:,2)-0.00;

% right foot
cop_vec=copr(:,1:2)-RFT_rc(:,1:2); % vector from com foot to cop
cop_dist=sqrt(sum(cop_vec.^2,2)); % distance of cop from com foot
cop_evec=[cop_vec(:,1)./cop_dist cop_vec(:,2)./cop_dist];
idx=cop_dist>foot_radius; % [m] Pythagoras]
copr(idx,1:2)=RFT_rc(idx,1:2) + foot_radius*cop_evec(idx,:);

% filter to make it a bit cleaner
%for iDim=1:2
%    copr(:,iDim)=filtfilt(B_filt,A_filt,copr(:,iDim));
%    copl(:,iDim)=filtfilt(B_filt,A_filt,copl(:,iDim));
%end

% total cop of both feet combined:
totMoment=cross(copr,grfr,2)+cross(copl,grfl,2);
grf_tot=grfl+grfr;
magGRF=dot(grf_tot,grf_tot,2);magGRF=repmat(magGRF,1,3);
r_coptot=cross(grf_tot,totMoment,2)./magGRF; % cop vector of right magnitude, perpendicular to force and moment

F_scale=(0-r_coptot(:,3))./(grf_tot(:,3)); % scale factor to find intersection of r_cop with ground plane
F_scale=repmat(F_scale,1,3);
r_coptot=r_coptot + F_scale.*grf_tot;

if doFig
    subplot(222);plot(copl(:,1),copl(:,2),'.'); hold on
    plot(copr(:,1),copr(:,2),'.');
    plot(r_coptot(:,1),r_coptot(:,2),'.');
    axis equal
    title('Corrected cop trajectories');legend('left','right','total')
    subplot(223);plot3(copl(:,1)-LFT_rc(:,1),copl(:,2)-LFT_rc(:,2),copl(:,3)-LFT_rc(:,3),'.'); hold on
    plot3(copr(:,1)-RFT_rc(:,1),copr(:,2)-RFT_rc(:,2),copr(:,3)-RFT_rc(:,3),'.'); title('cop wrt foot COM')
    axis equal
end
% calculate free moment
Mfree_l=[zeros(length(grfl),2) grml(:,3)-copl(:,1).*grfl(:,2)+copl(:,2).*grfl(:,1)];
Mfree_r=[zeros(length(grfl),2) grmr(:,3)-copr(:,1).*grfr(:,2)+copr(:,2).*grfr(:,1)];

if doFig
    fig=figure;
    fig.Name=['free moment check'];
    subplot(311);plot(grml(:,3),'r'); hold on; plot(grmr(:,3),'g');
    title('force plate vertical moment wrt mocap origin')
    subplot(312);plot(-copl(:,1).*grfl(:,2)+copl(:,2).*grfl(:,1),'r');
    hold on; plot(-copr(:,1).*grfr(:,2)+copr(:,2).*grfr(:,1),'g');
    title('moment of grf')
    subplot(313);plot(Mfree_l,'r');hold on;plot(Mfree_r,'g'); title('calculated free moments')
end
