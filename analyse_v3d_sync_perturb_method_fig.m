%% @reader:
% This code computes everything inverse dynamics, ie joint moments, forces
% and 6-DOF joint powers/work, COM power/work, segment mechanical
% energies, and residual forces and moments. All calculations are based on
% the following inputs from V3D: full segment kinematics, segment angular
% momenta (I*omega) and ground reaction forces. The GRF's are manipulated
% in this code to make them net zero in the horizontal, that is an optional
% step.
% ! All segment energies / power terms are expressed in a treadmill bound
% frame
% ! if you're lucky the dimension of the left hand side is behind each
% line, e.g. [m]: meters
% ! cntr-F @reader to find key points of explanation of some of the
% things in the code. Address questions to koenlemaire@posteo.net

clear all
%close all
clc

%% Switches, options and normalization constants:
do_residual_assignment=0; % assign residual forces other than the default
do_sanity_check=false; % produces some additional figs
animate=false; % make animation of raw data
% trial_nr_to_sanity_check=randperm(11,1); % choose from: [1..11]
trial_nr_to_sanity_check=12; % large1-2 med3-5 pref6-8 small9-11

% 'hard' makes mean(grf) hard zero in horizontal plane and makes grf and
% grm 0 during swing phase.
hard_force_correction=1;
make_non_dimensional=1; % non-dimensionalize some figures
do_force_filter=1; % filter forces

% warning in gait events detection function
give_event_detect_warning=1;

minColor=[0 0 1];
plusColor=[1 .3 0];

% delay control
delayVec=(-15:15)'; % samples delay in data
%delayVec=(-31:31)'; % samples delay in data


%% constants and parameters
belt_thickness=0.015; % [m]

% In this case the mocap reference frame has the same
% orientation but is positioned in front of the treadmill and a little
% to the left of the force plate coordinate system:
% mocap coordinate system:
% origin: forward left corner of left belt
%   wrt left belt origin; 1.626 m forward and 0.013 m left
% +X to right
% +Z up
% +Y forward
% vector of LEFT force plate origin in mocap reference frame;
r_left_FP_origin_in_mocap_frame=[0.013 -1.626 -belt_thickness];

% frequency at which forces and moments are filtered
force_filter_freq=20; % [Hz]

% non-dimensionalisation constants. approximate, only used to
% non-dimensionalize (should be individualized for multiple people)
L=0.94; % [m] leg length
g=9.81; % [m/s^2] gravitational acceleration
m=75; % [kg] mass [only used for normalization, doesn't matter]

if make_non_dimensional
    Pnorm = m*g*sqrt(g)*sqrt(L); % [W]
    Tnorm = m*g*L; % [Nm]
    Fnorm = m*g; % [N]
    angle_norm = 180/pi; % degree/rad
    time_norm=sqrt(L/g);
else
    Pnorm = 1; % [W]
    Tnorm = 1; % [Nm]
    Fnorm = 1; % [N]
    angle_norm = 1; % rad/degree
end

% allocation of belt speeds and trial types:
belt_speeds=[0.7 0.9 1.1 1.25 1.4 1.6 ... % const_sf and const_sf_2nd
    0.7 0.9 1.1 1.25 1.4 1.6 1.8 2.0 2.2 ... % pref
    0.9 1.1 0.7 0.9 1.1 0.7 0.9 1.1]; % large med small

% trial type pointers
idg_csf=1:6;
idg_pref=7:15;
airpump=16:23;
large=16:17;
med=18:20;
small=21:23;

% detection treshold for gait events (heelstrike and toe-off)
gait_detection_treshold=40; % [N]

% segment mass fractions. These data are taken from Visual 3D default
% settings. Fractions are for [foot shank thigh pelvis], where
% foot/shank/thigh is for ONE foot/shank/thigh
segment_mass_fractions=[0.0145 0.0465 0.1 0.142]'; % [kg]

%% path handling & loading data
% point d to the data location, and make sure the data directory is on the
% matlab path.

home_dir=mfilename('fullpath');
home_dir=home_dir(1:length(home_dir)-length(mfilename)); % full path of current file location
cd(home_dir) % move to directory containing file
addpath(genpath(cd))
%addpath(genpath('../airPump7_15/'))

%dataDir=dir('../soft_tissue_data/v3d_out_new');

for iFile=idg_pref(4)+2 % 1.25 m/s % 3:length(dataDir)
    % parse sanity check
    trial_nr=iFile-2;
    if do_sanity_check
        if ~all(trial_nr_to_sanity_check-(trial_nr))
            sanity_check=true;
        else
            sanity_check=false;
        end
    else
        sanity_check=false;
    end
    fname=invGrid_pref_04_synced.mat;%dataDir(iFile).name;
    mod=load(fname);
    
    %% Apply calibration and rotate force plate data to mocap frame:
    FP_raw=load([fname(1:end-7),'synced.mat']);
    grf_raw_V=FP_raw.forcePlateData_synced(:,[1:3,7:9]);
    grm_raw_V=FP_raw.forcePlateData_synced(:,[4:6,10:12]);
    [grfl_raw,grfr_raw,grml_raw,grmr_raw]=calibrateForces(grf_raw_V,grm_raw_V);
    
    %% parsing v3d variables
    % @ reader: V3D data includes kinematics (pos vel acc) of all segments
    % com and endpoints, angular kinematics and angular momenta of all
    % segments
    v3dNames=fields(mod);
    for i=1:length(v3dNames)
        mod.(v3dNames{i})=mod.(v3dNames{i}){:};% eval([v3dNames{i},'=v3dData.',v3dNames{i},'{1};'])
    end
    % we have learned that: the position data does not contain nan values
    % but the data requiring derivatives is missing datapoints at the
    % beginning and at the end (velocity and acceleration 3 and 4
    % datapoints respectively). We deal with this by chopping of the first
    % and last bit of the velocity data and then interpolating the force
    % data to these new values.
    
    % for all mocap data (INCLUDING THE TIME AXIS), but not the force
    % data we will chop off the first and last 5 samples to get rid of the
    % nans, then interpolate force data to new mocap time base
    n_mocap_frames=mod.frnr(end); % nr of mocap frames, note that after the next piece of code frnr will be altered!!
    
    for i=1:length(v3dNames)
        % logical indicating if this is a piece of data from the mocap
        % stuff, by cheking the length of data:
        if mod.(v3dNames{i})==n_mocap_frames% eval(['length(',v3dNames{i},')==n_mocap_frames;'])
            % now we chop a couple of frames:
            mod.(v3dNames{i})=mod.(v3dNames{i})(5:end-5,:); %eval([v3dNames{i},'=',v3dNames{i},'(5:end-5,:);'])
        else % we just do nothing
            continue
        end
    end
    
    % ! @reader for some reason
    mod.RTA_AngVel=(pi/180)*mod.RTA_AngVel;
    %% unify force and mocap time base
    interpMethod=['linear'];
    for iDim=1:3
        mod.grfl(:,iDim)=interp1(mod.frnr_analogTime,grfl_raw(:,iDim),mod.frnr_time,interpMethod); %[N]
        mod.grfr(:,iDim)=interp1(mod.frnr_analogTime,grfr_raw(:,iDim),mod.frnr_time,interpMethod); %[N]
        mod.grml(:,iDim)=interp1(mod.frnr_analogTime,grml_raw(:,iDim),mod.frnr_time,interpMethod); %[N]
        mod.grmr(:,iDim)=interp1(mod.frnr_analogTime,grmr_raw(:,iDim),mod.frnr_time,interpMethod); %[N]
    end
    
    mod.t=mod.frnr_time-mod.frnr_time(1); % [s] start at 0
    new_n_mocap_frames=length(mod.t);
    fs_mocap=FP_raw.debug.sample_rates.c3d;
    
    
    %% find heelstrike and to-off
    % @ reader here is a crude implementation of a tresholding method,
    % which appears to work well with this sample data, anything else
    % possible here
    
    % gait events in mocap time base.
    [hsl, tol, hsr, tor] = get_gait_events ...
        (mod.grfl(:,3),mod.grfr(:,3),gait_detection_treshold);
    
    % correct heelstrike indices; prepare to truncate data:
    switch trial_nr
        case idg_csf(1)
            hsl=hsl(1:end-3);
            hsr=hsr(1:end-3);
            tor=tor(1:end-3);
            tol=tol(1:end-3);
        case 7
            %tor=[tor(tor<9800); 9811; tor(tor>10000)];
        case 8
            %tor=[tor(tor<7800); 7996; tor(tor>8200)];
            %hsr=[hsr(hsr<7800); 8071; hsr(hsr>8200)];
        case 21
            hsr=[hsr(hsr<1500); 1790; hsr(hsr>1900)];
            tor=[tor(tor<1700); tor(tor>1900)];
            tor=[tor(tor<7100); tor(tor>7300)];
            tor=[tor(tor<11400); 11630; tor(tor>11800)];
            hsr=[hsr(hsr<11400); 11700; hsr(hsr>11800)];
        case 22
            hsl=hsl(2:end);
            tol=tol(2:end);
        case 23
            hsl=hsl(2:end);
            tol=tol(2:end);
    end
    
    % sanity check gait events
    [isbad] = check_gait_events (hsl, tol, hsr, tor, mod.grfl(:,3),mod.grfr(:,3), sanity_check);
    if isbad
        warning(['gait events in following trial nr are spurious: ',num2str(trial_nr)])
    end
    
    %% grf and grm correction based on swing phase
    % @reader this section substracts a line from the forces and moments
    % data. The line is constructed to pass through the forces and moments
    % during the periods the foot is off the ground
    [mod.grfl, mod.grfr, mod.grml, mod.grmr] = swing_phase_grf_correction...
        (mod.t,hsl,tol,hsr,tor,mod.grfl,mod.grfr,mod.grml,mod.grmr,sanity_check);
    % OBSOLETE: there was an idea to apply a rotation based on the mean force, in
    % such a way that the mean force vector is rotated to solely the
    % vertical direction. However, the rotation R*v1=v2 is not unique, ie
    % there is no unique R that satisfies this equation... Therefore, this
    % is a bad idea.
    % total grf:
    mod.grf_tot=mod.grfl+mod.grfr;
    %% transform force plate data to mocap reference frame
    [mod.grfl,mod.grfr,mod.grml,mod.grmr] = forcePlate_to_mocap_frame ...
        (mod.grfl,mod.grfr,mod.grml,mod.grmr,r_left_FP_origin_in_mocap_frame);
    
    %% filter forces
    if do_force_filter
        [B_filt,A_filt]=butter(4,force_filter_freq/FP_raw.debug.sample_rates.c3d/2);
        for iDim=1:3
            mod.grfl(:,iDim)=filtfilt(B_filt,A_filt,mod.grfl(:,iDim));
            mod.grfr(:,iDim)=filtfilt(B_filt,A_filt,mod.grfr(:,iDim));
            mod.grml(:,iDim)=filtfilt(B_filt,A_filt,mod.grml(:,iDim));
            mod.grmr(:,iDim)=filtfilt(B_filt,A_filt,mod.grmr(:,iDim));
        end
    end
    %% , truncate signals / perturbations to synchrony
    
        idx=hsl(1):hsl(end);
        hsl=hsl-idx(1)+1;
        tol=tol-idx(1)+1;
        hsr=hsr-idx(1)+1;
        tor=tor-idx(1)+1;
        oldMod=mod;
        old_idx=idx;
        
    for iDelay=1:length(delayVec)
        nr_samples_shift=delayVec(iDelay);
    
        % nr of strides
        n_strides=length(hsl);
        n_strides_all(trial_nr,1)=n_strides;
        
        
        % truncate signals and unify time:
        mod.t=oldMod.t(old_idx);
        mod.t=mod.t-mod.t(1); % start at 0
        mod.grfl=oldMod.grfl(old_idx,:);
        mod.grfr=oldMod.grfr(old_idx,:);
        %mod.copl_v3d=mod.copl_v3d(idx,:);
        %mod.copr_v3d=mod.copr_v3d(idx,:);
        mod.grml=oldMod.grml(old_idx,:);
        mod.grmr=oldMod.grmr(old_idx,:);
        %Mzl_raw_v3d=Mzl_raw_v3d(idx,:);
        %Mzr_raw_v3d=Mzr_raw_v3d(idx,:);
        
        %FP_raw.c3d_markerArray_synced=FP_raw.c3d_markerArray_synced(:,old_idx,:);
        
        % stride time
        stride_time(trial_nr,1)=mod.t(end)/n_strides;
        
        % step frequency
        step_freq(trial_nr,1)=2*n_strides/mod.t(end);
        
        % piece of code to truncate all signals such the we are left with and
        % integer number of strides in the dataset
        for i=1:length(v3dNames) % loop through all known variables coming out of v3d
            % truncate x if length(x) equals nr of mocap samples
            if length(oldMod.(v3dNames{i}))==new_n_mocap_frames % eval(['length(',char(v3dNames(i)),')==new_n_mocap_frames;'])
                mod.(v3dNames{i})=oldMod.(v3dNames{i})(old_idx+nr_samples_shift,:); % truncate
            else % do nothing
                continue
            end
        end
        % sanity check on grf:
        % conclusion: grf and grf_raw are exactly the same data (ie FP1 and
        % Force1 in LFT segment in v3d terms). heelstrike and toe-off detection
        % works OK. At this point data is unified in the time base
        
        %% grf correction to make swing phase grf=0 and horizontal mean(grf)=0
        if hard_force_correction
            [mod.grfl, mod.grfr, mod.grml, mod.grmr] = detrend_horizontal_grf...
                (mod.t,hsl,tol,hsr,tor,mod.grfl,mod.grfr,mod.grml,mod.grmr);
        end
        %% calculate COP and Free Moment from ground on subject
        % there are two methods for this
        % 1) principled method, given the GRF and GRM, the cop lies on a line
        % parallel to GRF at the appropriate distance from GRF. See calculation
        % of total COP ...
        % 2) use the formula provided in the book ...
        [mod.copl, mod.copr, mod.Mfree_l, mod.Mfree_r] = calculate_cop_Mfree ...
            (mod.grfl,mod.grfr,mod.grml,mod.grmr,mod.LFT_rc,mod.RFT_rc,sanity_check);
        
        %% animation of data
        if animate && sanity_check
            fig=figure;
            fig.Name=['animation'];
            markerArray=FP_raw.c3d_markerArray_synced;
            markerArray=markerArray/1000; % [m]
            markerArray(:,:,1)=-markerArray(:,:,1); % flip x direction
            markerArray(:,:,2)=-markerArray(:,:,2);
            for i=1:2:length(mod.grfl)
                % all markers:
                %tmp.XData=markerArray(:,i,1);
                %tmp.YData=markerArray(:,i,2);
                %tmp.ZData=markerArray(:,i,3);
                
                plot3(markerArray(:,i,1),markerArray(:,i,2),markerArray(:,i,3),'ko'); hold on
                
                plot3([mod.LSK_rc(i,1) mod.LSK_rc(i,1)],[mod.LSK_rc(i,2) mod.LSK_rc(i,2)],[mod.LSK_rc(i,3) mod.LSK_rc(i,3)],'ro')
                plot3([mod.RSK_rc(i,1) mod.RSK_rc(i,1)],[mod.RSK_rc(i,2) mod.RSK_rc(i,2)],[mod.RSK_rc(i,3) mod.RSK_rc(i,3)],'ro')
                
                plot3([mod.copl(i,1) mod.grfl(i,1)/Fnorm+mod.copl(i,1)],[mod.copl(i,2) mod.grfl(i,2)/Fnorm+mod.copl(i,2)],...
                    [mod.copl(i,3) mod.grfl(i,3)/Fnorm+mod.copl(i,3)],'r','linewidth',2)
                plot3(mod.copl(i,1),mod.copl(i,2),mod.copl(i,3),'ro')
                
                plot3([mod.copr(i,1) mod.grfr(i,1)/Fnorm+mod.copr(i,1)],[mod.copr(i,2) mod.grfr(i,2)/Fnorm+mod.copr(i,2) ],...
                    [mod.copr(i,3) mod.grfr(i,3)/Fnorm+mod.copr(i,3) ],'g','linewidth',2)
                tmp=plot3(mod.copr(i,1),mod.copr(i,2),mod.copr(i,3),'go');
                
                %xlabel x
                %ylabel y
                %zlabel z
                axis equal
                drawnow limitrate
                hold off
            end
        end
        %% COM power/work, subject mass
        % @reader: in this section COM power and work and related variables are
        % computed, and a grf correction is applied to make mean force in
        % horizontal plane 0. See function calc_com_power for details.
        
        % complete function call:
        % [mod.p_coml_belt, mod.p_comr_belt, v_com_belt, mod.v_com_lab, subj_mass, w_com_belt, dEkin_com_belt dEkin_com_lab] = calc_com_power(mod.grfl,mod.grfr,v_belt,t)
        [mod.p_coml, mod.p_comr, mod.v_com_belt, mod.a_com_lab, mod.v_com_lab, ...
            subj_mass(trial_nr,1), w_com(trial_nr,1), dEkin_com_belt(trial_nr,1), ...
            dmod.v_com_lab(trial_nr,:)] ...
            = calc_com_power_uncorrected(mod.grfl,mod.grfr,belt_speeds(trial_nr),mod.t);
        %dEkin_com_old(fcount,1)=.5*subj_mass(fcount)*dv*dv';% [J] kinetic energy change of COM
        mod.Pcom_tot=mod.p_coml+mod.p_comr;
        %W_com_net(trial_nr,1)=w_com(trial_nr,1); % [w] convert to average power
        %mod.grfl=mod.grfl_new;
        %mod.grfr=mod.grfr_new;
        
        %% assign segment masses
        % @reader mass is assigned here on a trial-by-trial basis, based on
        % mean(GRFZ) for that trial. Segment masses are adjusted accordingly
        % (but segment moments of intertia / angular momenta are not). Some
        % shenanigans with airpump weight, can skip that.
        
        % we do this here because at this point we have trial-to-trial subject
        % mass from our grf analysis (above).
        tot_mass=subj_mass(trial_nr,1);
        % tot_mass=74.77; % [kg] this is mean(subj_mass)
        if ~contains(fname,'invGrid') % we have a non pref, ie airpump trial
            mfoot=segment_mass_fractions(1)*(tot_mass-1)+0.5; % [kg] % correct for the airpump weight added
            mshank=segment_mass_fractions(2)*(tot_mass-1);
            mthigh=segment_mass_fractions(3)*(tot_mass-1);
            mpelvis=segment_mass_fractions(4)*(tot_mass-1);
            mRTA=(tot_mass)-2*(mfoot+mshank+mthigh)-mpelvis;
        else % no airpump trial
            mfoot=segment_mass_fractions(1)*tot_mass; % [kg]
            mshank=segment_mass_fractions(2)*tot_mass;
            mthigh=segment_mass_fractions(3)*tot_mass;
            mpelvis=segment_mass_fractions(4)*tot_mass;
            mRTA=tot_mass-2*(mfoot+mshank+mthigh)-mpelvis;
        end
        
        %% perturbations to angular momentum / linear momentum / synchrony
        % something I dream of in the future
        
        %% shadow inverse analysis: whole body linear and angular momentum balance
        % @reader the shadow inverse analysis is the critical part. What I call
        % 'shadow' is matlab based inverse dynamics calculation. In principle
        % this can be compared to V3D based calculation; I have done that in
        % the past and that checked out. Currently this comparison has non-zero
        % difference because I meddled with the GRF's / masses (see above), but
        % you can trust this to be accurate.
        % In this section I first compute linear momentum balance for the body
        % as a whole, which yields the residual force. Then compute angular
        % momentum balance about the COM, which yields the residual moment. I
        % will use these later to apply them on different segments (see
        % explanation in next section). This calculation is optional (but
        % cool!)
        
        % linear momentum balance equation: sum of forces equals rate of change
        % of linear momentum of the COM
        
        % sum of forces, left hand side (lhs)
        mod.grf_tot=mod.grfl+mod.grfr;
        com_F_grav=zeros(size(mod.RTA_rcdd)); % gravitational force
        com_F_grav(:,3)=-g*tot_mass;
        mod.Fsum_lhs=mod.grf_tot+com_F_grav; % [N]
        
        % rate of change of linear momentum COM (right hand side, rhs)
        mod.pdot_rhs = mfoot*(mod.RFT_rcdd + mod.LFT_rcdd) + ... % feet
            mshank*(mod.RSK_rcdd + mod.LSK_rcdd) + ... % shanks
            mthigh*(mod.RTH_rcdd + mod.LTH_rcdd) + ... % thighs
            mpelvis*mod.RPV_rcdd + ... % pelvis
            mRTA*mod.RTA_rcdd; % [kg*m/s2] trunk
        
        
        % residual: difference between right and left hand side
        mod.Fres_body=mod.pdot_rhs - mod.Fsum_lhs; % [N]
        mod.Fres_mean(trial_nr,:)=trapz(mod.t,mod.Fres_body)/(mod.t(end)-mod.t(1)); % [N]
        mod.Fres_mag=sqrt(sum(mod.Fres_body.^2,2));
        mod.Fres_mean_mag(trial_nr,1)=trapz(mod.t,mod.Fres_mag)/(mod.t(end)-mod.t(1)); % [N]
        % Angular momentum balance about COM: sum of moments about COM equals
        % the rate of change of angular momentum about COM
        
        % position of COM
        mod.r_COM = (mfoot*(mod.RFT_rc + mod.LFT_rc) + ... % feet
            mshank*(mod.RSK_rc + mod.LSK_rc) + ... % shanks
            mthigh*(mod.RTH_rc + mod.LTH_rc) + ... % thighs
            mpelvis*mod.RPV_rc + ... % pelvis
            mRTA*mod.RTA_rc)/tot_mass; % trunk
        
        % IMPROVE here using just measured GRM
        % sum of moments about COM: moments from GRF (r x F) plus 'free' moment
        % from ground on body (about vertical axis) (left hand side of angmom
        % balance equation):
        Mfree_tot=mod.Mfree_l+mod.Mfree_r;
        Mgrfl_wrt_COM=cross(mod.copl-mod.r_COM,mod.grfl,2);
        Mgrfr_wrt_COM=cross(mod.copr-mod.r_COM,mod.grfr,2);
        Mgrf_tot_wrt_COM=Mgrfr_wrt_COM+Mgrfl_wrt_COM;
        Msum_wrt_COM_lhs2=Mfree_tot+Mgrf_tot_wrt_COM; % [Nm]
        
        mod.Msum_wrt_COM_lhs=mod.grml+mod.grmr+...
            cross(-mod.r_COM,mod.grfl)+cross(-mod.r_COM,mod.grfr); % [Nm]
        % Angular momentum of one segment about body COM equals angular
        % momentum about segment COM ('spin' angular momentum) + the moment of
        % the linear momentum of that segment wrt to the body COM ('orbital'
        % angular momentum). Total angular momentum equals the sum of the
        % segment angular momenta about the body COM.
        
        % sum of segment's spin angular momenta:
        H_spin=mod.LFT_AngMom+mod.RFT_AngMom+mod.LSK_AngMom+mod.RSK_AngMom+...
            mod.LTH_AngMom+mod.RTH_AngMom+mod.RPV_AngMom+mod.RTA_AngMom; % [kg*m2/s] rotational part for each segment
        % sum of 'orbital' angular momentum about body COM for each segment (ie
        % cross product of position vector wrt body COM and linear momentum
        % vector of segment COM)
        H_orbital=cross(mod.RFT_rc-mod.r_COM,mfoot*mod.RFT_rcd,2)+...
            cross(mod.LFT_rc-mod.r_COM,mfoot*mod.LFT_rcd,2)+...
            cross(mod.RSK_rc-mod.r_COM,mshank*mod.RSK_rcd,2)+...
            cross(mod.LSK_rc-mod.r_COM,mshank*mod.LSK_rcd,2)+...
            cross(mod.RTH_rc-mod.r_COM,mthigh*mod.RTH_rcd,2)+...
            cross(mod.LTH_rc-mod.r_COM,mthigh*mod.LTH_rcd,2)+...
            cross(mod.RPV_rc-mod.r_COM,mpelvis*mod.RPV_rcd,2)+...
            cross(mod.RTA_rc-mod.r_COM,mRTA*mod.RTA_rcd,2); % [kg*m2/s]
        % total angular momentum about COM, sum of spin and orbital angular
        % momentum for each segment
        mod.H_seg_wrt_COM=H_spin+H_orbital;
        
        % rate of change angular momentum about COM; the time derivative of
        % total body angular momentum (right hand side). (@XFU you do something
        % else for this I believe, which is numerically more precise. I've
        % found however that this method results in near-numerical precision
        % error wrt to V3D based calculations)
        mod.Hdot_wrt_COM_rhs=(diff1d(mod.H_seg_wrt_COM',fs_mocap))'; %[kg*m2/s2]
        
        % residual: error in angular momentum balance equation
        mod.Mres_body_wrt_COM=mod.Hdot_wrt_COM_rhs-mod.Msum_wrt_COM_lhs; %[Nm]
        mod.Mres_mag=sqrt(sum(mod.Mres_body_wrt_COM.^2,2));
        if sanity_check
            %keyboard
            %mean(Mres2_wrt_COM)
            fig=figure;
            fig.Name=['angular momentum balance'];
            subplot(321);plot(mod.t,mod.Msum_wrt_COM_lhs)
            title('Sum of moments wrt COM')
            xlabel('Time [s]')
            ylabel('Moment [Nm]')
            subplot(323);plot(mod.t,mod.Hdot_wrt_COM_rhs)
            title('rate of change of angular momentum wrt COM')
            xlabel('Time [s]')
            ylabel('dH/dt [kg*m^2/s^2')
            subplot(325);plot(mod.t,mod.Mres_body_wrt_COM)
            title('Residual moment about COM')
            xlabel('Time [s]')
            ylabel('Moment [Nm]')
            
            subplot(322);plot(mod.t,mod.Fsum_lhs)
            title('Rate of change of linear momentum')
            xlabel('Time [s]')
            ylabel('rate of linear momentum [kg*m/s^2]')
            subplot(324);plot(mod.t,mod.pdot_rhs)
            title('Sum of forces')
            xlabel('Time [s]')
            ylabel('force [N]')
            subplot(326);plot(mod.t,mod.Fres_body)
            title('Residual force')
            xlabel('Time [s]')
            ylabel('force [N]')
            
            fig=figure;
            fig.Name=['angular momentum contributions about vertical'];
            subplot(411);plot(mod.t,[Mgrfl_wrt_COM(:,3) Mgrfr_wrt_COM(:,3) Mgrf_tot_wrt_COM(:,3)])
            title('moment of GRF wrt COM')
            xlabel('Time [s]')
            ylabel('Moment [Nm]')
            legend('left','right','total')
            subplot(412);plot(mod.t,[mod.Mfree_l(:,3) mod.Mfree_r(:,3) Mfree_tot(:,3)])
            title('Free moment, force couple')
            xlabel('Time [s]')
            ylabel('Moment [Nm]')
            legend('left','right','total')
            subplot(413);plot(mod.t,mod.Hdot_wrt_COM_rhs(:,3))
            title('rate of change of angular momentum wrt COM')
            xlabel('Time [s]')
            ylabel('dH/dt [kg*m^2/s^2')
            subplot(414);plot(mod.t,mod.Mres_body_wrt_COM(:,3))
            title('Residual moment about COM')
            xlabel('Time [s]')
            ylabel('Moment [Nm]')
            
        end
        %% assignment of residual force and moment
        % @reader !!this is tricky!!, but it turns out that you can assign the
        % residual force on any segment of the body, on any place of that
        % segment, provided that you choose your residuals in such a way that
        % the total linear and angular momentum balance equations (see above)
        % are satisfied. To me this is a crucial point, because the joint
        % moments and powers that you end up with depend on this choice!! When
        % you start at the bottom (conventional inverse dynamics) you
        % implicitly assign all the residuals at the end of the chain, on the
        % head. Alternatively you could start at the top and assign these
        % residuals to the feet, or actually distribute them in any way over
        % all segments! Again, how they are distributed influences the outcome
        % of joint moment and power calculations. I actually think very few
        % people in the biomechanics community realize this .... And this was
        % exactly my point of contention with Art ...
        
        % for now we'll assign residuals proportional to body mass of each
        % segment. If we make the (partial) residual force act each segment's
        % COM, than the net residual force act on the body COM, meaning that we
        % can next assign the residual moments at will, while still satisfying
        % linear and angular momentum balance.
        if do_residual_assignment
            mod.Fres_foot = mfoot/tot_mass*mod.Fres_body;
            mod.Fres_shank = mshank/tot_mass*mod.Fres_body;
            mod.Fres_thigh = mthigh/tot_mass*mod.Fres_body;
            mod.Fres_pelvis = mpelvis/tot_mass*mod.Fres_body;
            mod.Fres_trunk = mRTA/tot_mass*mod.Fres_body;
            
            mod.Mres_foot = mfoot/tot_mass*mod.Mres_body_wrt_COM;
            mod.Mres_shank = mshank/tot_mass*mod.Mres_body_wrt_COM;
            mod.Mres_thigh = mthigh/tot_mass*mod.Mres_body_wrt_COM;
            mod.Mres_pelvis = mpelvis/tot_mass*mod.Mres_body_wrt_COM;
            mod.Mres_trunk = mRTA/tot_mass*mod.Mres_body_wrt_COM;
        else
            mod.Fres_foot=zeros(size(mod.grfl));
            mod.Fres_shank=zeros(size(mod.grfl));
            mod.Fres_thigh=zeros(size(mod.grfl));
            mod.Fres_pelvis=zeros(size(mod.grfl));
            mod.Fres_trunk=zeros(size(mod.grfl));
            
            mod.Mres_foot=zeros(size(mod.grfl));
            mod.Mres_shank=zeros(size(mod.grfl));
            mod.Mres_thigh=zeros(size(mod.grfl));
            mod.Mres_pelvis=zeros(size(mod.grfl));
            mod.Mres_trunk=zeros(size(mod.grfl));
        end
        % potentially we can now play around with distributing the residual
        % moment in a different way, and show that this affects the results!
        %mod.Mres_thigh = .5*mRTA/tot_mass*Mres2_wrt_COM;
        %mod.Mres_trunk = mpelvis/tot_mass*Mres2_wrt_COM;
        %mod.Mres_pelvis = mthigh/tot_mass*Mres2_wrt_COM;
        %% shadow inverse analysis: dynamics
        % @reader this is the meat of the inverse dynamics calculation,
        % calculating all joint moments and force by doing angular and linear
        % momentum balance for each segment, sequentially.
        % this section is thoroughly checked (KKL 9/2020); meaning the shadow
        % calculated joint forces and torques yielded the same values (within
        % numerical precision) as those provided by v3d, if the subject mass
        % and grf's are unaltered from v3d output and all residuals where
        % assigned to the final segment, ie if you just ignore/omit the
        % residuals in the bottom equations you end up with 'regular' inverse
        % dynamics analysis. At each step we calculate the 'proximal' force /
        % torque, which is believed to act at the proximal point of the segment
        % (r_prox, coming from V3D), and which is based on the 'downstream'
        % forces and moments, and the kinematics of the segment.
        
        % --forces. We will do this the dumb way, ie working from the ground up
        % rather than solving a system of equations simultaneously.
        
        
        % forces of gravity on each segment:
        F_grav_FT=zeros(size(mod.LFT_rcdd)); % [N] foot
        F_grav_FT(:,3)=-g*mfoot;
        
        F_grav_SK=zeros(size(mod.LFT_rcdd)); % shank
        F_grav_SK(:,3)=-g*mshank;
        
        F_grav_TH=zeros(size(mod.LFT_rcdd)); % etc ...
        F_grav_TH(:,3)=-g*mthigh;
        
        F_grav_RPV=zeros(size(mod.LFT_rcdd));
        F_grav_RPV(:,3)=-g*mpelvis;
        
        F_grav_RTA=zeros(size(mod.LFT_rcdd));
        F_grav_RTA(:,3)=-g*mRTA;
        
        % feet linear momentum balance
        mod.LFT_F_prox=mfoot*mod.LFT_rcdd - (F_grav_FT + mod.grfl + mod.Fres_foot); % [N] brackets needed here
        mod.RFT_F_prox=mfoot*mod.RFT_rcdd - (F_grav_FT + mod.grfr + mod.Fres_foot);
        
        % shanks
        mod.LSK_F_prox=mshank*mod.LSK_rcdd - F_grav_SK + mod.LFT_F_prox - mod.Fres_shank; % no brackets because LFT_F_prox has minus sign on left-hand side
        mod.RSK_F_prox=mshank*mod.RSK_rcdd - F_grav_SK + mod.RFT_F_prox - mod.Fres_shank;
        
        % thighs
        mod.LTH_F_prox=mthigh*mod.LTH_rcdd - F_grav_TH + mod.LSK_F_prox - mod.Fres_thigh;
        mod.RTH_F_prox=mthigh*mod.RTH_rcdd - F_grav_TH + mod.RSK_F_prox - mod.Fres_thigh;
        
        % Hip
        mod.RPV_F_prox=mpelvis*mod.RPV_rcdd - F_grav_RPV + mod.LTH_F_prox + mod.RTH_F_prox - mod.Fres_pelvis;
        
        % Trunk.
        % @reader Fres on this segment should be zero because we distributed the
        % total Fres over all segments! This is a validation of our
        % calculations. If you omit Fres from all segments than the residual
        % calculated on the trunk should (and will) equal the residual computed
        % from the whole body linear momentum balance above
        mod.Fres_on_trunk_COM=mRTA*mod.RTA_rcdd - F_grav_RTA + mod.RPV_F_prox - mod.Fres_trunk; % Newton 2nd law
        
        % --- moments. again doing this the dumb way ...
        % angular momentum balance: sum of torques about segment com must
        % equal rate of change of angular momentum of segment:
        % feet (see comments for foot, pattern repeats moving up the body)
        % again one can omit the residuals in a regular calculation
        
        % Torque of force plate about foot com (this will be replaced by torque of
        % lower lying segment force when moving up)
        T_FP_LFT=mod.grml+cross(-mod.LFT_rc,mod.grfl); % [Nm]
        T_FP_RFT=mod.grmr+cross(-mod.RFT_rc,mod.grfr); % [Nm]
        
        % old version based on COP
        %T_F_dist_LFT = cross(mod.copl-mod.LFT_rc,mod.grfl,2); % [Nm]
        %T_F_dist_RFT = cross(mod.copr-mod.RFT_rc,mod.grfr,2);
        
        % torque of force from shank on foot (computed above) about com foot:
        T_F_prox_LFT = cross(mod.LFT_r_prox-mod.LFT_rc,mod.LFT_F_prox,2); % [Nm]
        T_F_prox_RFT = cross(mod.RFT_r_prox-mod.RFT_rc,mod.RFT_F_prox,2);
        
        if sanity_check
            figure;
            subplot(321);plot(mod.grfl)
            subplot(322);plot(mod.grfr)
            subplot(323);plot(T_FP_LFT)
            subplot(324);plot(T_FP_RFT)
            subplot(325);plot(T_F_prox_LFT)
            subplot(326);plot(T_F_prox_RFT)
        end
        
        % torque from shank on foot, based on angular momentum balance for this
        % segment. This is the result of one iteration and is used moving up
        % the chain.
        mod.LFT_T_prox=(diff1d(mod.LFT_AngMom',fs_mocap))' - (T_FP_LFT+T_F_prox_LFT+mod.Mres_foot); % [Nm]
        mod.RFT_T_prox=(diff1d(mod.RFT_AngMom',fs_mocap))' - (T_FP_RFT+T_F_prox_RFT+mod.Mres_foot);
        
        % old version using Mfree and moment of GRF, based on COP, this is
        % verified to yield the same result!!
        % mod.LFT_T_prox2=(diff1d(mod.LFT_AngMom',fs_mocap))' - (T_F_dist_LFT+T_F_prox_LFT+mod.Mfree_l+mod.Mres_foot); % [Nm]
        % mod.RFT_T_prox2=(diff1d(mod.RFT_AngMom',fs_mocap))' - (T_F_dist_RFT+T_F_prox_RFT+mod.Mfree_r+mod.Mres_foot);
        
        % shanks
        T_F_dist_LSK = cross(mod.LFT_r_prox-mod.LSK_rc,-mod.LFT_F_prox,2);
        T_F_dist_RSK = cross(mod.RFT_r_prox-mod.RSK_rc,-mod.RFT_F_prox,2);
        
        T_F_prox_LSK = cross(mod.LSK_r_prox-mod.LSK_rc,mod.LSK_F_prox,2);
        T_F_prox_RSK = cross(mod.RSK_r_prox-mod.RSK_rc,mod.RSK_F_prox,2);
        
        mod.LSK_T_prox=(diff1d(mod.LSK_AngMom',fs_mocap))' - (T_F_dist_LSK+T_F_prox_LSK-mod.LFT_T_prox+mod.Mres_shank);
        mod.RSK_T_prox=(diff1d(mod.RSK_AngMom',fs_mocap))' - (T_F_dist_RSK+T_F_prox_RSK-mod.RFT_T_prox+mod.Mres_shank);
        
        % thighs
        T_F_dist_LTH = cross(mod.LSK_r_prox-mod.LTH_rc,-mod.LSK_F_prox,2);
        T_F_dist_RTH = cross(mod.RSK_r_prox-mod.RTH_rc,-mod.RSK_F_prox,2);
        
        T_F_prox_LTH = cross(mod.LTH_r_prox-mod.LTH_rc,mod.LTH_F_prox,2);
        T_F_prox_RTH = cross(mod.RTH_r_prox-mod.RTH_rc,mod.RTH_F_prox,2);
        
        mod.LTH_T_prox=(diff1d(mod.LTH_AngMom',fs_mocap))' - (T_F_dist_LTH+T_F_prox_LTH-mod.LSK_T_prox+mod.Mres_thigh);
        mod.RTH_T_prox=(diff1d(mod.RTH_AngMom',fs_mocap))' - (T_F_dist_RTH+T_F_prox_RTH-mod.RSK_T_prox+mod.Mres_thigh);
        
        % pelvis (dist is towards feet, prox is towards head)
        % @reader here the forces and torques from both legs acting at the hips
        % need to be added, this requires some bookkeeping but is not in
        % principle different from any of the above.
        T_F_dist_LPV = cross(mod.LTH_r_prox-mod.RPV_rc,-mod.LTH_F_prox,2);
        T_F_dist_RPV = cross(mod.RTH_r_prox-mod.RPV_rc,-mod.RTH_F_prox,2);
        
        T_F_prox_PV = cross(mod.RPV_r_prox-mod.RPV_rc,mod.RPV_F_prox,2);
        
        mod.RPV_T_prox=(diff1d(mod.RPV_AngMom',fs_mocap))' - ...
            (T_F_dist_LPV+T_F_dist_RPV+T_F_prox_PV -mod.LTH_T_prox -mod.RTH_T_prox +mod.Mres_pelvis);
        
        % trunk (dist towards legs, prox (=res) towards head)
        % @reader here it gets tricky again: if the total Mres has been
        % distributed already Mres1 should be zero. In the traditional approach
        % (all residuals ignored thus far) the residual moment on the final
        % segment (Mres1) equals the Mres from whole body angular momentum
        % balance (see above) PLUS the moment of the residual force about the
        % COM. This took me a while to realize ..... but it turns out to be
        % true. If ignored, the residual force in the force equations acts at
        % the COM of the trunk, and thus has a moment about the COM of the
        % body, which appears in the net residual moment acting on the trunk..
        % You see that all this gets complicated very quickly ... ask @Koen
        % Lemaire for detailed explanation
        
        T_F_dist_RTA = cross(mod.RPV_r_prox-mod.RTA_rc,-mod.RPV_F_prox,2); % acting from pelvis on trunk
        % angular momentum balance for the trunk
        mod.Mres_about_trunk_COM=(diff1d(mod.RTA_AngMom',fs_mocap))' -T_F_dist_RTA +mod.RPV_T_prox -mod.Mres_trunk;
        % residual moment about body COM
        mod.Mres_trunk_about_body_COM = mod.Mres_about_trunk_COM - cross(mod.r_COM-mod.RTA_rc,mod.Fres_on_trunk_COM,2);
        
        % sanity check: shadow force and moment calculations. If the residuals
        % are all distributed on the segments the final residuals on the trunk
        % equal zero, otherwise they will reflect the total residuals in the
        % whole body angular/linear momentum balance.
        if sanity_check % these should correspond to v3d output if nothing is perturbed
            fig=figure;
            fig.Name=['sanity check of force and moment balance'];
            subplot(221);plot(mod.t,[mod.Fres_on_trunk_COM])
            ylabel('residual force [N]')
            title('Fres on trunk, after applying forces')
            subplot(222);plot(mod.t,[mod.Fres_body])
            title('Fres using total rate of change of linear momentum and grf')
            
            subplot(223); plot(mod.t,mod.Mres_about_trunk_COM); hold on
            title('Mres on trunk, after applying moments')
            ylabel('Residual moment')
            xlabel('Time [s]')
            
            subplot(224); plot(mod.t,mod.Mres_body_wrt_COM); hold on
            title('Mres using total angular momentum and grf/grm')
            ylabel('Residual moment')
            xlabel('Time [s]')
        end
        
        %% shadow inverse analysis: (rate of change of) energies
        % @reader here the segment mechanical energies are calculated based on
        % V3D kinematic variables, as the sum of potential gravitational and
        % kinetic energy. See notes throughout.
        % !!! ALLL !!! power/work/energy calculations are made in
        % treadmill frame.
        % !! Note that the v3d provided segment energies / total body energy is
        % in lab frame, and can thus not be used.
        % !! Note that soft tissue power depends on reference frame, and that
        % the treadmill frame IS the reference frame for which this variable
        % assumes its 'true' value, in the sense that this is the frame in
        % which the soft tissue power corresponds to heat change associated
        % with it.
        
        % define belt speed vector:
        v_belt=zeros(size(mod.LFT_rcd));
        v_belt(:,2)=belt_speeds(trial_nr);
        
        % total mechanical energy of each segment: Ekin + Epot. Note the
        % addition of vbelt in the calculation of the translational kinetic
        % energy. Note that the change in total mechanical energy corresponds
        % to the 6 DOF powers exerted on each segment (if these were to be
        % calculated). Note that the potential energy term is positive for
        % positve vertical coordinates, as it should be.
        mod.Emech_LFT=0.5*(mfoot*sum((mod.LFT_rcd+v_belt).^2,2) + sum(mod.LFT_AngMom.*mod.LFT_AngVel,2)) + mfoot*mod.LFT_rc(:,3)*g; % g is a positive number here!!
        mod.Emech_RFT=0.5*(mfoot*sum((mod.RFT_rcd+v_belt).^2,2) + sum(mod.RFT_AngMom.*mod.RFT_AngVel,2)) + mfoot*mod.RFT_rc(:,3)*g; % [J]
        
        mod.Emech_LSK=0.5*(mshank*sum((mod.LSK_rcd+v_belt).^2,2) + sum(mod.LSK_AngMom.*mod.LSK_AngVel,2)) + mshank*mod.LSK_rc(:,3)*g;
        mod.Emech_RSK=0.5*(mshank*sum((mod.RSK_rcd+v_belt).^2,2) + sum(mod.RSK_AngMom.*mod.RSK_AngVel,2)) + mshank*mod.RSK_rc(:,3)*g;
        
        mod.Emech_LTH=0.5*(mthigh*sum((mod.LTH_rcd+v_belt).^2,2) + sum(mod.LTH_AngMom.*mod.LTH_AngVel,2)) + mthigh*mod.LTH_rc(:,3)*g;
        mod.Emech_RTH=0.5*(mthigh*sum((mod.RTH_rcd+v_belt).^2,2) + sum(mod.RTH_AngMom.*mod.RTH_AngVel,2)) + mthigh*mod.RTH_rc(:,3)*g;
        
        mod.Emech_RPV=0.5*(mpelvis*sum((mod.RPV_rcd+v_belt).^2,2) + sum(mod.RPV_AngMom.*mod.RPV_AngVel,2)) + mpelvis*mod.RPV_rc(:,3)*g;
        
        mod.Emech_RTA=0.5*(mRTA*sum((mod.RTA_rcd+v_belt).^2,2) + sum(mod.RTA_AngMom.*mod.RTA_AngVel,2)) + mRTA*mod.RTA_rc(:,3)*g;
        
        % model total mechanical energy in treadmill frame
        mod.Emech_model=mod.Emech_LFT+mod.Emech_RFT+mod.Emech_LSK+mod.Emech_RSK+mod.Emech_LTH+mod.Emech_RTH+mod.Emech_RPV+mod.Emech_RTA; % [J]
        
        % model rate of energy change:
        mod.Emech_modeld=(diff1d(mod.Emech_model',fs_mocap))';
        mod.delta_Emech_model(trial_nr,1)=(mod.Emech_model(end)-mod.Emech_model(1));
        
        % Whole body peripheral energy (independent of reference frame, because
        % relative to com bound reference frame)
        mod.Eper_rot=0.5*(sum([ ... % rotational peripheral energy [J]
            mod.LFT_AngMom.*mod.LFT_AngVel mod.RFT_AngMom.*mod.RFT_AngVel ...
            mod.LSK_AngMom.*mod.LSK_AngVel mod.RSK_AngMom.*mod.RSK_AngVel ...
            mod.LTH_AngMom.*mod.LTH_AngVel mod.RTH_AngMom.*mod.RTH_AngVel ...
            mod.RPV_AngMom.*mod.RPV_AngVel mod.RTA_AngMom.*mod.RTA_AngVel ],2));
        
        mod.Eper_trans = ... % translational peripheral energy [J]
            0.5*mfoot*sum([((mod.RFT_rcd-mod.v_com_lab).^2) (mod.LFT_rcd-mod.v_com_lab).^2],2)+ ...
            0.5*mshank*sum([((mod.RSK_rcd-mod.v_com_lab).^2) (mod.LSK_rcd-mod.v_com_lab).^2],2)+ ...
            0.5*mthigh*sum([((mod.RTH_rcd-mod.v_com_lab).^2) (mod.LTH_rcd-mod.v_com_lab).^2],2)+ ...
            0.5*mpelvis*sum([(mod.RPV_rcd-mod.v_com_lab).^2],2) + ...
            0.5*mRTA*sum([(mod.RTA_rcd-mod.v_com_lab).^2],2);
        
        %clear mod.v_com_lab
        
        % total peripheral energy, ie sum of kinetic energy of all mass
        % relative to body COM
        mod.Eper=mod.Eper_rot+mod.Eper_trans;
        
        mod.Pper=(diff1d(mod.Eper',fs_mocap))';
        mod.Eper_net(trial_nr,1)=(mod.Eper(end)-mod.Eper(1));
        
        mod.P_com_tot=mod.p_coml+mod.p_comr;
        mod.P_com_PLUS_per=mod.P_com_tot+mod.Pper;
        
        idx=mod.P_com_PLUS_per>0;
        W_comPlusPer(trial_nr,:)=[trapz(mod.t(idx),mod.P_com_PLUS_per(idx)) trapz(mod.t(~idx),mod.P_com_PLUS_per(~idx))]/n_strides;
        
        W_com_PLUS_per(trial_nr,1)=trapz(mod.t,mod.P_com_PLUS_per); % [J]
        
        %W_com_PLUS_per(trial_nr,1)=(mod.Eper(end)-mod.Eper(1)+w_com(trial_nr,1))/n_strides;
        if sanity_check
            fig=figure;
            fig.Name=['comparison of Pcom+per to Emechdot'];
            subplot(211)
            plot(mod.t,mod.p_coml,'r',mod.t(tol),mod.p_coml(tol),'r+',mod.t(hsl),mod.p_coml(hsl),'ro'); hold on
            plot(mod.t,mod.p_comr,'b',mod.t(tor),mod.p_comr(tor),'b+',mod.t(hsr),mod.p_comr(hsr),'bo')
            title('left (red) and right (blue) com power')
            xlabel('time [s]')
            ylabel('Pcom [W]')
            subplot(212);plot(mod.t,[mod.P_com_PLUS_per mod.Emech_modeld])
            legend('Pcom + per','rate of change of sum of seg energies')
            xlabel('time [s]')
            ylabel('Power [W]')
        end
        %% shadow inverse analysis: joint power/work, residual powers
        % @reader here all 6-DOF joint powers are calculated, in such a way the
        % the 6-DOF joint powers, plus powers of residual force and moment (if
        % these are assigned) are equal to the rate of change of the segment
        % energy. The calculation is relatively straightfoward, see 6-DOF power
        % function for details.
        % !!! ALLL !!! power/work/energy calculations are made in
        % treadmill frame. See comments at segment energies above.
        % ! Note that distal foot power represents the 6-DOF power between the
        % ground and the foot; in some sense, one could argue that this should
        % be part of 'joint' power, in the sense that it should be needed to
        % satisfy the energy balance equation of the body
        % ! the 'residual power' is the difference between the sum of joint
        % powers and the total rate of change in whole body mechanical energy.
        % If the latter is computed from the sum of segment energies I (KKL)
        % would call this residual power. If the total body rate of change is
        % calculated as COM + PER power, we call this soft tissue power.
        
        %feet:
        mod.LFT_rcd_belt=v_belt+mod.LFT_rcd;
        mod.RFT_rcd_belt=v_belt+mod.RFT_rcd;
        
        mod.LFT_vcop=cross(mod.LFT_AngVel,mod.copl-mod.LFT_rc);
        mod.RFT_vcop=cross(mod.RFT_AngVel,mod.copr-mod.RFT_rc);
        
        mod.r_cop_LFT=mod.copl-mod.LFT_rc;
        mod.r_cop_RFT=mod.copr-mod.RFT_rc;
        
        mod.v_copl_loc=v_belt+mod.LFT_rcd + cross(mod.LFT_AngVel,mod.copl-mod.LFT_rc);
        mod.v_copr_loc=v_belt+mod.RFT_rcd + cross(mod.RFT_AngVel,mod.copr-mod.RFT_rc);
        
        v_origin_l_local=v_belt+mod.LFT_rcd + cross(mod.LFT_AngVel,-mod.LFT_rc);
        v_origin_r_local=v_belt+mod.RFT_rcd + cross(mod.RFT_AngVel,-mod.RFT_rc);
        
        % split out for debugging purposes ...
        mod.Pgrfl=sum(mod.v_copl_loc.*mod.grfl,2);
        mod.Pgrfr=sum(mod.v_copr_loc.*mod.grfr,2);
        
        mod.Pgrfl_vec=mod.v_copl_loc.*mod.grfl;
        mod.Pgrfr_vec=mod.v_copr_loc.*mod.grfr;
        
        mod.Pgrml=sum(mod.LFT_AngVel.*mod.Mfree_l,2);
        mod.Pgrmr=sum(mod.RFT_AngVel.*mod.Mfree_r,2);
        
        % old version, completely equivalent!!
        mod.Pdist_lfoot2=mod.Pgrfl + mod.Pgrml;
        mod.Pdist_rfoot2=mod.Pgrfr + mod.Pgrmr;
        
        mod.Pdist_lfoot=dot(v_origin_l_local,mod.grfl,2) + dot(mod.LFT_AngVel,mod.grml,2);
        mod.Pdist_rfoot=dot(v_origin_r_local,mod.grfr,2) + dot(mod.RFT_AngVel,mod.grmr,2);
        
        if sanity_check
            fig=figure;fig.Name=['distal foot power contributions. trial_nr=',num2str(trial_nr)];
            subplot(321);plot(mod.v_copl_loc(:,2),'r'); hold on;plot(mod.v_copr_loc(:,2),'g');
            title('cop velocity wrt foot COM')
            subplot(322);plot(sum(mod.v_copl_loc.*mod.grfl,2),'r'); hold on; plot(sum(mod.v_copr_loc.*mod.grfr,2),'g')
            title('power from grf')
            subplot(323);plot(sum(mod.LFT_AngVel.*mod.Mfree_l,2),'r'); hold on; plot(sum(mod.RFT_AngVel.*mod.Mfree_r,2),'g')
            title('power from grm')
            subplot(324);plot(mod.Pdist_lfoot,'r'); hold on; plot(mod.Pdist_rfoot,'g')
            title('total power')
            
            % power balance for feet:
            mod.Emech_LFTd=diff1d(mod.Emech_LFT',fs_mocap)';
            mod.Emech_RFTd=diff1d(mod.Emech_RFT',fs_mocap)';
            
            % velocity of ankle joint on shank body (in treadmill frame!):
            % power of left foot on left shank:
            P_LSK_on_LFT = sum((mod.LFT_rd_prox+v_belt).*mod.LFT_F_prox,2) + sum(mod.LFT_AngVel.*mod.LFT_T_prox,2);
            P_RSK_on_RFT = sum((mod.RFT_rd_prox+v_belt).*mod.RFT_F_prox,2) + sum(mod.RFT_AngVel.*mod.RFT_T_prox,2);
            
            % power of residual force and moment on left shank (in treadmill
            % frame)
            Pres_on_LFT=sum((mod.LFT_rcd+v_belt).*mod.Fres_foot,2) + sum(mod.LFT_AngVel.*mod.Mres_foot,2);
            Pres_on_RFT=sum((mod.RFT_rcd+v_belt).*mod.Fres_foot,2) + sum(mod.RFT_AngVel.*mod.Mres_foot,2);
            
            % total power on left shank
            Ptot_LFT=mod.Pdist_lfoot+P_LSK_on_LFT+Pres_on_LFT;
            Ptot_RFT=mod.Pdist_rfoot+P_RSK_on_RFT+Pres_on_RFT;
            
            subplot(325);plot(Ptot_LFT); hold on; plot(mod.Emech_LFTd,'--')
            title('power balance left')
            legend('total power','dEmechdt')
            subplot(326);plot(Ptot_RFT); hold on; plot(mod.Emech_RFTd,'--')
            title('power balance right')
            legend('total power','dEmechdt')
        end
        
        
        % Ankle, Knee, Hip, left side:
        % note that separating joint power in a '3-DOF' and '6-DOF' part makes
        % limited sense for the case that the mechanics are 6-DOF ...
        % see Calc_6DOFJointPower for details
        mod.LeftAnklePower= Calc_6DOFJointPower(mod.LSK_AngVel,mod.LFT_AngVel,...
            mod.LSK_rcd,mod.LSK_rc,mod.LFT_rd_prox,mod.LFT_r_prox,mod.LFT_F_prox,mod.LFT_T_prox); % [W]
        mod.LeftKneePower = Calc_6DOFJointPower(mod.LTH_AngVel,mod.LSK_AngVel,...
            mod.LTH_rcd,mod.LTH_rc,mod.LSK_rd_prox,mod.LSK_r_prox,mod.LSK_F_prox,mod.LSK_T_prox);
        mod.LeftHipPower = Calc_6DOFJointPower(mod.RPV_AngVel,mod.LTH_AngVel,...
            mod.RPV_rcd,mod.RPV_rc,mod.LTH_rd_prox,mod.LTH_r_prox,mod.LTH_F_prox,mod.LTH_T_prox);
        
        % right side
        mod.RightAnklePower = Calc_6DOFJointPower(mod.RSK_AngVel,mod.RFT_AngVel,...
            mod.RSK_rcd,mod.RSK_rc,mod.RFT_rd_prox,mod.RFT_r_prox,mod.RFT_F_prox,mod.RFT_T_prox); %[W]
        mod.RightKneePower = Calc_6DOFJointPower(mod.RTH_AngVel,mod.RSK_AngVel,...
            mod.RTH_rcd,mod.RTH_rc,mod.RSK_rd_prox,mod.RSK_r_prox,mod.RSK_F_prox,mod.RSK_T_prox);
        mod.RightHipPower = Calc_6DOFJointPower(mod.RPV_AngVel,mod.RTH_AngVel,...
            mod.RPV_rcd,mod.RPV_rc,mod.RTH_rd_prox,mod.RTH_r_prox,mod.RTH_F_prox,mod.RTH_T_prox);
        
        
        % pelvis-trunk
        mod.S5L1Power=Calc_6DOFJointPower(mod.RTA_AngVel,mod.RPV_AngVel,...
            mod.RTA_rcd,mod.RTA_rc,mod.RPV_rd_prox,mod.RPV_r_prox,mod.RPV_F_prox,mod.RPV_T_prox); % [W]
        
        % residual:
        % for the trunk, the sum of all powers must equal the rate of change of
        % its kinetic energy: Pjoint + Pres = d_Emech_dt
        
        
        % velocity of mod.RPV_r_prox (S5L1 joint) on trunk body:
        v_S5L1_RTA=v_belt + mod.RTA_rcd + cross(mod.RTA_AngVel,mod.RPV_r_prox-mod.RTA_rc,2); % [m/s]
        % power of pelvis on trunk:
        mod.P_pelvis_trunk=sum(v_S5L1_RTA.*-mod.RPV_F_prox,2) + sum(mod.RTA_AngVel.*-mod.RPV_T_prox,2);
        Pres_trunk_energy_balance=diff1d(mod.Emech_RTA',fs_mocap)' - mod.P_pelvis_trunk; % residual power based on energy balance
        Pres_trunk_direct=sum(mod.RTA_AngVel.*mod.Mres_about_trunk_COM,2) + sum([(v_belt+mod.RTA_rcd).*mod.Fres_on_trunk_COM],2); % residual power based on Tres and Fres
        
        % we have learned that:
        % -RTA_F_prox just satisfies Newton for trunk segment, meaning de-facto
        % that -mod.RPV_F_prox does not act on the trunk segment. Instead,
        % RTA_F_prox acts on the trunk segment, and part of that might be in
        % mod.RPV_F_prox, or something.
        % -The proximal end of RTA is the part closest to the pelvis, the
        % distal end is at the head.
        % -Indeed, Newton for whole body and for only trunk segment yield the
        % same Fres
        
        
        %% sum of joint powers and avg joint powers + power balance sanity check:
        % @reader here we just sum power / work values calculated above
        % sum of joint power and joint work
        % positive / negative joint work calculations
        idxl=mod.LeftAnklePower>0;
        idxr=mod.RightAnklePower>0;
        W_ankle(iDelay,:)=[(trapz(mod.t(idxl),mod.LeftAnklePower(idxl))+trapz(mod.t(idxr),mod.RightAnklePower(idxr)))/2 ...
             (trapz(mod.t(~idxl),mod.LeftAnklePower(~idxl))+trapz(mod.t(~idxr),mod.RightAnklePower(~idxr)))/2];
        
        idxl=mod.LeftKneePower>0;
        idxr=mod.RightKneePower>0;
        W_knee(iDelay,:)=[(trapz(mod.t(idxl),mod.LeftKneePower(idxl))+trapz(mod.t(idxr),mod.RightKneePower(idxr)))/2 ...
             (trapz(mod.t(~idxr),mod.RightKneePower(~idxr))+trapz(mod.t(~idxl),mod.LeftKneePower(~idxl)))/2];
        
        idxl=mod.LeftHipPower>0;
        idxr=mod.RightHipPower>0;
        W_hip(iDelay,:)=[(trapz(mod.t(idxl),mod.LeftHipPower(idxl))+trapz(mod.t(idxr),mod.RightHipPower(idxr)))/2 ...
            (trapz(mod.t(~idxr),mod.RightHipPower(~idxr))+trapz(mod.t(~idxl),mod.LeftHipPower(~idxl)))/2];
        
        idxl=mod.Pdist_lfoot>0;
        idxr=mod.Pdist_rfoot>0;
        W_fdist(trial_nr,:)=[trapz(mod.t(idxl),mod.Pdist_lfoot(idxl)) trapz(mod.t(~idxl),mod.Pdist_lfoot(~idxl))...
            trapz(mod.t(idxr),mod.Pdist_rfoot(idxr)) trapz(mod.t(~idxr),mod.Pdist_rfoot(~idxr))]/n_strides;
        
        Wres_net(trial_nr,1)=trapz(mod.t,Pres_trunk_energy_balance);
        
        mod.sumJointPower_no_feet=mod.RightAnklePower+mod.RightKneePower+mod.RightHipPower+ ...
            mod.LeftAnklePower+mod.LeftKneePower+mod.LeftHipPower+mod.S5L1Power; % [W]
        
        idx=mod.sumJointPower_no_feet>0;
        W_sumJoint_no_feet(trial_nr,:)=[trapz(mod.t(idx),mod.sumJointPower_no_feet(idx)) trapz(mod.t(~idx),mod.sumJointPower_no_feet(~idx))]/n_strides;
        
        mod.sumJointPower_with_feet=mod.sumJointPower_no_feet+mod.Pdist_lfoot+mod.Pdist_rfoot;
        
        AnkleWork(iDelay,1)=trapz(mod.t,mod.LeftAnklePower+mod.RightAnklePower);
        KneeWork(iDelay,1)=trapz(mod.t,mod.LeftKneePower+mod.RightKneePower);
        HipWork(iDelay,1)=trapz(mod.t,mod.LeftHipPower+mod.RightHipPower); %[J]
        S5L1Work(trial_nr,1)=trapz(mod.t,mod.S5L1Power);
        W_dist_foot(trial_nr,1)=trapz(mod.t,mod.Pdist_lfoot+mod.Pdist_rfoot);
        
        AnkleWork(iDelay,1)
        
        sumJointWork_no_feet(trial_nr,1)=trapz(mod.t,mod.sumJointPower_no_feet);
        sumJointWork_with_feet(trial_nr,1)=trapz(mod.t,mod.sumJointPower_with_feet);
        
        if sanity_check
            fig=figure;
            fig.Name=['sanity check on power calculations'];
            LFT_LSK_AngVel=mod.LFT_AngVel-mod.LSK_AngVel;
            left_ankle_rot_power_check=LFT_LSK_AngVel.*mod.LFT_T_prox*pi/180; % left foot angular velocity expressed wrt shank reference frame DOT proximal torque on left foot
            subplot(221)
            plot(sum(left_ankle_rot_power_check,2))
            title('left ankle rotational power')
            legend('v3d in LSK coord','v3d in LAB coord','v3d as scalar','from joint omega dot Torque')
            ylabel('power [W]')
            
            % first sanity check roational kinetic energy calculation, erot
            % is equal to .5 * angular momentum innner prodoct angular velocity:
            left_shank_erot2=.5*sum(mod.LSK_AngMom.*mod.LSK_AngVel,2);
            % sanity check power calculations for left shank:
            % the sum of all power should equal the rate of change of kinetic
            % energy of the segment
            % segment rate of change of mechanical energy:
            
            mod.Emech_LSKd=(diff1d(mod.Emech_LSK',fs_mocap))';
            % segment angular velocities are in rad/s according to:
            % https://c-motion.com/v3dwiki/index.php?title=KINETIC_KINEMATIC:_AngVel
            
            % segment energies are in Joules, as suggested by:
            % https://c-motion.com/v3dwiki/index.php?title=TRANSLATIONAL_ENERGY_SCALAR
            % powers:
            
            % velocity of ankle joint on shank body (in treadmill frame!):
            v_ankle_LSK=mod.LSK_rcd+v_belt + cross(mod.LSK_AngVel,mod.LFT_r_prox-mod.LSK_rc,2); % [m/s]
            % power of left foot on left shank:
            P_lfoot_lshank = sum(v_ankle_LSK.*-mod.LFT_F_prox,2) + sum(mod.LSK_AngVel.*-mod.LFT_T_prox,2);
            
            % power of left thigh on left shank (in treadmill frame!):
            P_lthigh_lshank=sum((mod.LSK_rd_prox+v_belt).*mod.LSK_F_prox,2) + sum(mod.LSK_AngVel.*mod.LSK_T_prox,2);
            
            % power of residual force and moment on left shank (in treadmill
            % frame)
            P_res_lshank=sum((mod.LSK_rcd+v_belt).*mod.Fres_shank,2) + sum(mod.LSK_AngVel.*mod.Mres_shank,2);
            
            % total power on left shank
            Ptot_lshank=P_lfoot_lshank + P_lthigh_lshank + P_res_lshank;
            
            subplot(222);plot(mod.t,[mod.left_shank_erot left_shank_erot2])
            title('v3d based erot vs 0.5*angmom*omega')
            legend('v3d based erot','0.5*angmom*omega')
            ylabel('Kinetic Energy [J]')
            
            subplot(223);plot(mod.t,[mod.Emech_LSKd Ptot_lshank])
            title('power balance left shank')
            legend('rate of change Emech','sum of powers')
            xlabel('time [s]')
            ylabel('power [W]')
            
            subplot(224);plot(mod.t,[mod.Pdist_lfoot mod.Pdist_rfoot])
            title('distal foot power')
            legend('left','right')
            xlabel('time [s]')
            ylabel('power [W]')
        end
        
        
        %% soft tissue power
        % @reader this is the soft tissue power and work per Art's calculation,
        % pretty straightforward
        mod.Psoft_with_feet=mod.P_com_PLUS_per-mod.sumJointPower_with_feet; % [W]
        mod.Psoft_no_feet=mod.P_com_PLUS_per-mod.sumJointPower_no_feet;
        idx=mod.Psoft_no_feet>0;
        W_soft(trial_nr,:)=[trapz(mod.t(idx),mod.Psoft_no_feet(idx)) trapz(mod.t(~idx),mod.Psoft_no_feet(~idx))]/n_strides;
        
        Wsoft_with_feet(trial_nr,1)=trapz(mod.t,mod.Psoft_with_feet); % [J]
        Wsoft_no_feet(trial_nr,1)=trapz(mod.t,mod.Psoft_no_feet);
        %% collision work
        for iStride=1:length(hsl)-1
            idx=(hsl(iStride):hsl(iStride+1));
            tmp=mod.p_coml(idx);
            zci=find(diff(sign(tmp)));
            if tmp(1)<0 % we already start at collision!
                idx_col=1:zci(1);
            else
                idx_col=zci(1):zci(2);
            end
            coll_work_l_stride(iStride,1)=trapz(tmp(idx_col))/FP_raw.debug.sample_rates.c3d;
        end
        coll_work_l(trial_nr,1)=mean(coll_work_l_stride);
        
        for iStride=1:length(hsr)-1
            idx=(hsr(iStride):hsr(iStride+1));
            tmp=mod.p_comr(idx);
            zci=find(diff(sign(tmp)));
            if tmp(1)<0 % we already start at collision!
                idx_col=1:zci(1);
            else
                idx_col=zci(1):zci(2);
            end
            coll_work_r_stride(iStride,1)=trapz(tmp(idx_col))/FP_raw.debug.sample_rates.c3d;
        end
        coll_work_r(trial_nr,1)=mean(coll_work_r_stride);
        
        clearvars *_stride
        %% average over strides wrt left heelstrike
        base=linspace(0,1,101); % 101 base units
        modNames=fieldnames(mod);
        for iVar=1:length(modNames)
            current_var=mod.(modNames{iVar}); % current variable that we are dealing with
            [m,n]=size(current_var);
            if m~=length(mod.grfl)
                continue % we don't even want to look at this
            end
            stride_nr=0; % we do this to be able to skip strides in a convenient way
            for iStride=1:length(hsl)-1
                stride_nr=stride_nr+1;
                idx=(hsl(iStride):hsl(iStride+1));
                norm_base=linspace(0,1,length(idx));
                for iDim=1:n % looping through all variables
                    current_var_stride(stride_nr,iDim,:)=interp1(norm_base,current_var(idx,iDim),base);
                end
                if n>1 % vectorial variable
                    avgMod.(modNames{iVar})(:,:,iDelay)=squeeze(mean(current_var_stride,1))';
                else % scalar variable
                    avgMod.(modNames{iVar})(:,iDelay)=squeeze(mean(current_var_stride,1))';
                end
                
            end
            clearvars current_var_stride
        end
        
        %% individual trial figures
        % @reader the rest is making figures
        if sanity_check
            % com power
            % 3x3 raw data left
            fig=figure;
            fig.Name=['left normal 3x3 ',mod.FILE_NAME];
            subplot(3,3,1)
            plot(mod.t,mod.Left_Ankle_Angle)
            title('ankle angle')
            xlabel('Time [s]')
            ylabel('Angle [deg]')
            legend('saggital','frontal','transverse')
            subplot(3,3,2)
            plot(mod.t,mod.Left_Knee_Angle)
            title('knee angle')
            xlabel('Time [s]')
            ylabel('Angle [deg]')
            subplot(3,3,3)
            plot(mod.t,mod.Left_Hip_Angle)
            title('hip angle')
            xlabel('Time [s]')
            ylabel('Angle [deg]')
            
            subplot(3,3,4)
            plot(mod.t,mod.LFT_T_prox)
            title('ankle moment')
            xlabel('Time [s]')
            ylabel('moment [Nm]')
            subplot(3,3,5)
            plot(mod.t,mod.LSK_T_prox)
            xlabel('Time [s]')
            ylabel('moment [Nm]')
            title('knee moment')
            subplot(3,3,6)
            plot(mod.t,mod.LTH_T_prox)
            xlabel('Time [s]')
            ylabel('moment [Nm]')
            title('hip moment')
            
            subplot(3,3,7)
            plot(mod.t,[mod.LeftAnklePower])
            xlabel('Time [s]')
            ylabel('Power [W]')
            title('ankle power')
            
            subplot(3,3,8)
            plot(mod.t,[mod.LeftKneePower])
            xlabel('Time [s]')
            ylabel('Power [W]')
            title('knee power')
            subplot(3,3,9)
            plot(mod.t,[mod.LeftHipPower])
            xlabel('Time [s]')
            ylabel('Power [W]')
            title('hip power')
            
            % right leg
            fig=figure;
            fig.Name=['right normal 3x3 ',mod.FILE_NAME];
            subplot(3,3,1)
            plot(mod.t,mod.Right_Ankle_Angle)
            title('ankle angle')
            xlabel('Time [s]')
            ylabel('Angle [deg]')
            legend('saggital','frontal','transverse')
            subplot(3,3,2)
            plot(mod.t,mod.Right_Knee_Angle)
            title('knee angle')
            xlabel('Time [s]')
            ylabel('Angle [deg]')
            subplot(3,3,3)
            plot(mod.t,mod.Right_Hip_Angle)
            title('hip angle')
            xlabel('Time [s]')
            ylabel('Angle [deg]')
            
            subplot(3,3,4)
            plot(mod.t,mod.RFT_T_prox)
            title('ankle moment')
            xlabel('Time [s]')
            ylabel('moment [Nm]')
            subplot(3,3,5)
            plot(mod.t,mod.RSK_T_prox)
            xlabel('Time [s]')
            ylabel('moment [Nm]')
            title('knee moment')
            subplot(3,3,6)
            plot(mod.t,mod.RTH_T_prox)
            xlabel('Time [s]')
            ylabel('moment [Nm]')
            title('hip moment')
            
            subplot(3,3,7)
            plot(mod.t,[mod.RightAnklePower])
            xlabel('Time [s]')
            ylabel('Power [W]')
            title('ankle power')
            
            subplot(3,3,8)
            plot(mod.t,[mod.RightKneePower])
            xlabel('Time [s]')
            ylabel('Power [W]')
            title('knee power')
            subplot(3,3,9)
            plot(mod.t,[mod.RightHipPower])
            xlabel('Time [s]')
            ylabel('Power [W]')
            title('hip power')
        end
    end % loop over files!!!
end % loop over delays!!!

%% PAPERFIG -- color definitions
base=linspace(0,100,101);
Lcolor=[0.8 0 0];
Rcolor=[0 .49 .38];
Ncolor=[0 0 0];
alphaMod=linspace(0.4,1,3);alphaMod=fliplr(alphaMod);
TimeLabel=['Time (\% gait cycle)'];
lblx=0.05;
lbly=1.05;


%% PAPERFIG 3x3
tlim=[-.15 .25];
alim=[-1.5 .5];
plim=[-.14 .25];
flim=[-.1 1.3];
figDelays=[1,ceil(length(delayVec)/2),length(delayVec)];
DIM=1;
for shiftCon=1
    fig=figure(23);
    axlimit=[0 1.6 -400 400];
    switch shiftCon
        case 1
            shiftnr=50;
            shift_name='circshift';
        case 2
            shiftnr=0;
            shift_name=[];
    end
    fig.Name=['3x3 ',shift_name];
    for iDelay=1:3
        leftCol=[Lcolor alphaMod(iDelay)];
        rightCol=[Rcolor alphaMod(iDelay)];
        switch iDelay
            case 1
                totCol=[plusColor alphaMod(iDelay)];
            case 2
                totCol=[Ncolor alphaMod(iDelay)];        
            case 3
                totCol=[minColor alphaMod(iDelay)];
        end
        idx_tmp=figDelays(iDelay);
        %base=base_stride(idx_tmp,:);
        
        % angle
        subplot(3,4,2);
        
        %plot(base,-avgMod.Left_Ankle_Angle(:,DIM,idx_tmp)/angle_norm,'color',leftCol,'linestyle','-.'); hold on
        %plot(base,circshift(-avgMod.Right_Ankle_Angle(:,DIM,idx_tmp)/angle_norm,shiftnr),'color',rightCol); hold on
        if iDelay==2
        plot(base,(-avgMod.Left_Ankle_Angle(:,DIM,idx_tmp)+circshift(-avgMod.Right_Ankle_Angle(:,DIM,idx_tmp),shiftnr))/2/angle_norm,'-','color',Ncolor); hold on
        
        %plot(18+[0 bar_length(idx_tmp)],-0.23*[4-iSpeed 4-iSpeed]-.6,'color',totCol)
        %plot(18+[0 bar_length(idx_tmp)],-0.23*[4-iSpeed 4-iSpeed]-.6,'color',totCol)
        %text(18+bar_length(idx_tmp)+8-1.7*iSpeed,-0.23*(4-iSpeed)-.6,speed_labels(iSpeed),'fontsize',8)
        %if iSpeed==3
        %text(10+bar_length(idx_tmp)/2,-0.23*(4-iSpeed)-.4,'0.2 s','fontsize',8)
        %end       
        title('Ankle')
        ylabel('Angle (rad)')
        ylim(alim)
        box off
        
        subplot(3,4,3)
        %plot(base,avgMod.Left_Knee_Angle(:,DIM,idx_tmp)/angle_norm,'color',leftCol,'linestyle','-.'); hold on
        %plot(base,circshift(avgMod.Right_Knee_Angle(:,DIM,idx_tmp)/angle_norm,shiftnr),'color',rightCol); hold on
        plot(base,(circshift(avgMod.Right_Knee_Angle(:,DIM,idx_tmp),shiftnr)+avgMod.Left_Knee_Angle(:,DIM,idx_tmp))/2/angle_norm,'-','color',Ncolor); hold on
        
        title('Knee')
        ylim(alim)
        box off
        
        tmp=subplot(3,4,4);
        yyaxis left
        tmp.YColor=[0 0 0];
        %plot(base,-avgMod.Left_Hip_Angle(:,DIM,idx_tmp)/angle_norm,'color',leftCol,'linestyle','-.'); hold on
        %plot(base,circshift(-avgMod.Right_Hip_Angle(:,DIM,idx_tmp)/angle_norm,shiftnr),'color',rightCol); hold on
        plot(base,(-avgMod.Left_Hip_Angle(:,DIM,idx_tmp)+circshift(-avgMod.Right_Hip_Angle(:,DIM,idx_tmp),shiftnr))/2/angle_norm,'-','color',Ncolor); hold on
        
        title('Hip')
        ylim(alim)
        tmp2=tmp.YLim;
        yyaxis right
        ylabel('(deg)')
        tmp.YColor=[0 0 0];
        tmp.YLim=tmp2*angle_norm;
        box off
        end 
        % moment
        subplot(3,4,6)
        %plot(base,-avgMod.LFT_T_prox(:,DIM,idx_tmp)/Tnorm,'color',leftCol,'linestyle','-.');  hold on
        %plot(base,circshift(-avgMod.RFT_T_prox(:,DIM,idx_tmp)/Tnorm,shiftnr),'color',rightCol);  hold on
        plot(base,(circshift(-avgMod.RFT_T_prox(:,DIM,idx_tmp),shiftnr)-avgMod.LFT_T_prox(:,DIM,idx_tmp))/2/Tnorm,'color',totCol);  hold on
        if iDelay==1
            annotation('textarrow',[0.20 0.20],[0.45 0.56],'String','$+ \Delta t \,$')
        end
        ylabel('Moment ($MgL$)')
        ylim(tlim)
        box off
        
        subplot(3,4,7)
        %plot(base,avgMod.LSK_T_prox(:,DIM,idx_tmp)/Tnorm,'color',leftCol,'linestyle','-.'); hold on
        %plot(base,circshift(avgMod.RSK_T_prox(:,DIM,idx_tmp)/Tnorm,shiftnr),'color',rightCol); hold on
        plot(base,(circshift(avgMod.RSK_T_prox(:,DIM,idx_tmp),shiftnr)+avgMod.LSK_T_prox(:,DIM,idx_tmp))/2/Tnorm,'color',totCol); hold on
        
        ylim(tlim)
        box off
        
        tmp=subplot(3,4,8);
        yyaxis left
        tmp.YColor=[0 0 0];
        %plot(base,-avgMod.LTH_T_prox(:,DIM,idx_tmp)/Tnorm,'color',leftCol,'linestyle','-.'); hold on
        %plot(base,circshift(-avgMod.RTH_T_prox(:,DIM,idx_tmp)/Tnorm,shiftnr),'color',rightCol); hold on
        plot(base,(circshift(-avgMod.RTH_T_prox(:,DIM,idx_tmp),shiftnr)-avgMod.LTH_T_prox(:,DIM,idx_tmp))/2/Tnorm,'-','color',totCol); hold on
        
        ylim(tlim)
        tmp2=tmp.YLim;
        yyaxis right
        ylabel('(N m)')
        tmp.YColor=[0 0 0];
        tmp.YLim=tmp2*Tnorm;
        box off
        
        % power
        subplot(3,4,10)
        %plot(base,[avgMod.LeftAnklePower(:,idx_tmp)]/Pnorm,'color',leftCol,'linestyle','-.'); hold on
        %plot(base,circshift([avgMod.RightAnklePower(:,idx_tmp)]/Pnorm,shiftnr),'color',rightCol); hold on
        plot(base,(circshift(avgMod.RightAnklePower(:,idx_tmp),shiftnr)+avgMod.LeftAnklePower(:,idx_tmp))/2/Pnorm,'color',totCol); hold on
        ylabel('Power ($Mg^{1.5} L^{0.5}$)')
        ylim(plim)
        xlabel(TimeLabel)
        box off
        
        subplot(3,4,11)
        %plot(base,[avgMod.LeftKneePower(:,idx_tmp)]/Pnorm,'color',leftCol,'linestyle','-.'); hold on
        %plot(base,circshift([avgMod.RightKneePower(:,idx_tmp)]/Pnorm,shiftnr),'color',rightCol); hold on
        plot(base,(circshift(avgMod.RightKneePower(:,idx_tmp),shiftnr)+avgMod.LeftKneePower(:,idx_tmp))/2/Pnorm,'color',totCol); hold on
        
        ylim(plim)
        xlabel(TimeLabel)
        box off
        
        tmp=subplot(3,4,12);
        yyaxis left
        tmp.YColor=[0 0 0];
        %plot(base,[avgMod.LeftHipPower(:,idx_tmp)]/Pnorm,'color',leftCol,'linestyle','-.'); hold on
        %plot(base,circshift([avgMod.RightHipPower(:,idx_tmp)]/Pnorm,shiftnr),'color',rightCol); hold on
        plot(base,(circshift(avgMod.RightHipPower(:,idx_tmp),shiftnr)+[avgMod.LeftHipPower(:,idx_tmp)])/2/Pnorm,'-','color',totCol); hold on
        
        ylim(plim)
        xlabel(TimeLabel)
        tmp2=tmp.YLim;
        yyaxis right
        ylabel('(W)')
        tmp.YColor=[0 0 0];
        tmp.YLim=tmp2*Pnorm;
        box off
    end
end

%% PAPERFIG residual forces (time)
figure(69)

for iDelay=1:3
    switch iDelay
        case 1
            totCol=[plusColor alphaMod(iDelay)];
        case 2
            totCol=[Ncolor alphaMod(iDelay)];
        case 3
            totCol=[minColor alphaMod(iDelay)];
    end
        idx_tmp=figDelays(iDelay);
        %base=base_stride(idx_tmp,:);
        
        % FRES
        tmp=subplot(2,4,3);
        yyaxis left
        tmp.YColor=[0 0 0];
        plot(base,avgMod.Fres_mag(:,idx_tmp)/Fnorm,'-','color',totCol); hold on
        if iDelay==1
        %text(tmp,lblx,lbly,'E','units','normalized')
        title('Walking data')
        end
        box off
        ylabel('Residual force ($Mg$)')
        tmp2=tmp.YLim;
        xlabel(TimeLabel)
        
        yyaxis right
        tmp.YColor=[0 0 0];
        tmp.YLim=tmp2*Fnorm;
        box off
        ylabel('(N)')
        
        % FRES
        tmp=subplot(2,4,4);
        yyaxis left
        tmp.YColor=[0 0 0];
        plot(base,avgMod.Mres_mag(:,idx_tmp)/Tnorm,'-','color',totCol); hold on
        if iDelay==1
        %text(tmp,lblx,lbly,'F','units','normalized')
        end
        
        %plot(18+[0 bar_length(idx_tmp)],-0.23*[4-iSpeed 4-iSpeed]-.6,'color',totCol)
        %plot(18+[0 bar_length(idx_tmp)],-0.23*[4-iSpeed 4-iSpeed]-.6,'color',totCol)
        %text(18+bar_length(idx_tmp)+8-1.7*iSpeed,-0.23*(4-iSpeed)-.6,speed_labels(iSpeed),'fontsize',8)
        %if iSpeed==3
        %text(10+bar_length(idx_tmp)/2,-0.23*(4-iSpeed)-.4,'0.2 s','fontsize',8)
        %end       
        box off
        ylabel('Residual moment ($Mg$)')
        tmp2=tmp.YLim;
        xlabel(TimeLabel)
        
        yyaxis right
        tmp.YColor=[0 0 0];
        tmp.YLim=tmp2*Tnorm;
        box off
        ylabel('($\mathrm{N} \, \mathrm{m}$)')
        
end
%% PAPERFIG -- RMSE and WORK
zero_delay_idx=find(delayVec==0);

figure(50)
for iDelay=1:length(delayVec)
    % angles
    tmp=squeeze(avgMod.Left_Ankle_Angle(:,:,iDelay)-avgMod.Left_Ankle_Angle(:,:,zero_delay_idx));
    ankle_angle_error(iDelay)=sqrt(mean(sum(tmp.^2,2)));
    tmp=squeeze(avgMod.Left_Knee_Angle(:,:,iDelay)-avgMod.Left_Knee_Angle(:,:,zero_delay_idx));
    knee_angle_error(iDelay)=sqrt(mean(sum(tmp.^2,2)));
    tmp=squeeze(avgMod.Left_Hip_Angle(:,:,iDelay)-avgMod.Left_Hip_Angle(:,:,zero_delay_idx));
    hip_angle_error(iDelay)=sqrt(mean(sum(tmp.^2,2)));
    % torques
    tmp1=squeeze(avgMod.LFT_T_prox(:,:,iDelay)-avgMod.LFT_T_prox(:,:,zero_delay_idx));
    tmp2=squeeze(avgMod.RFT_T_prox(:,:,iDelay)-avgMod.RFT_T_prox(:,:,zero_delay_idx));
    LFT_T_prox_error(iDelay)=(sqrt(mean(sum(tmp1.^2,2))) + sqrt(mean(sum(tmp2.^2,2))))/2;
    tmp=squeeze(avgMod.LTH_T_prox(:,:,iDelay)-avgMod.LTH_T_prox(:,:,zero_delay_idx));
    LTH_T_prox_error(iDelay)=sqrt(mean(sum(tmp.^2,2)));
    tmp=squeeze(avgMod.LSK_T_prox(:,:,iDelay)-avgMod.LSK_T_prox(:,:,zero_delay_idx));
    LSK_T_prox_error(iDelay)=sqrt(mean(sum(tmp.^2,2)));
    % powers
    tmp1=squeeze(avgMod.LeftAnklePower(:,iDelay)-avgMod.LeftAnklePower(:,zero_delay_idx));
    tmp2=squeeze(avgMod.RightAnklePower(:,iDelay)-avgMod.RightAnklePower(:,zero_delay_idx));
    LeftAnklePower_error(iDelay)=sqrt(mean(tmp1.^2));
    AnklePowerRMSE(iDelay)=(sqrt(mean(tmp1.^2)) + sqrt(mean(tmp2.^2)))/2;
    tmp=squeeze(avgMod.LeftKneePower(:,iDelay)-avgMod.LeftKneePower(:,zero_delay_idx));
    LeftKneePower_error(iDelay)=sqrt(mean(tmp.^2));
    tmp=squeeze(avgMod.LeftHipPower(:,iDelay)-avgMod.LeftHipPower(:,zero_delay_idx));
    LeftHipPower_error(iDelay)=sqrt(mean(tmp.^2));
    
    mean_Fres(iDelay)=sqrt(mean(avgMod.Fres_mag(:,iDelay).^2))/Fnorm;
    mean_Mres(iDelay)=sqrt(mean(avgMod.Mres_mag(:,iDelay).^2))/Tnorm;    
end
LFT_T_prox_error_rel=LFT_T_prox_error/sqrt(mean(sum(squeeze(avgMod.LFT_T_prox(:,:,zero_delay_idx)).^2,2)));
LTH_T_prox_error_rel=LTH_T_prox_error/sqrt(mean(sum(squeeze(avgMod.LTH_T_prox(:,:,zero_delay_idx)).^2,2)));
LSK_T_prox_error_rel=LSK_T_prox_error/sqrt(mean(sum(squeeze(avgMod.LSK_T_prox(:,:,zero_delay_idx)).^2,2)));

AnklePowerRMSE_rel=AnklePowerRMSE/mean([ sqrt(mean(squeeze(avgMod.LeftAnklePower(:,zero_delay_idx)).^2)) sqrt(mean(squeeze(avgMod.RightAnklePower(:,zero_delay_idx)).^2)) ]);

HipPowerRMSE_rel=LeftHipPower_error/sqrt(mean(squeeze(avgMod.LeftHipPower(:,zero_delay_idx)).^2));

exp_rsme_data.delays=delayVec/960;
exp_rsme_data.moment=[LFT_T_prox_error];
exp_rsme_data.Tnorm=Tnorm;
exp_rsme_data.power=[AnklePowerRMSE];
exp_rsme_data.Pnorm=Pnorm;
exp_rsme_data.mean_Fres=mean_Fres*Fnorm;
exp_rsme_data.Fnorm=Fnorm;

save('exp_rsme_data','exp_rsme_data')

tmp=subplot(326);
%yyaxis left
tmp.YColor=[0 0 0];
plot(delayVec,W_ankle(:,1)/n_strides/Tnorm,'k-'); hold on
plot(delayVec,W_knee(:,1)/n_strides/Tnorm,'k--'); hold on
plot(delayVec,W_hip(:,1)/n_strides/Tnorm,'k:','linewidth',1.1); hold on

plot(delayVec,W_ankle(:,2)/n_strides/Tnorm,'k-'); hold on
plot(delayVec,W_knee(:,2)/n_strides/Tnorm,'k--'); hold on
plot(delayVec,W_hip(:,2)/n_strides/Tnorm,'k:','linewidth',1.1); hold on

xlabel('Delay $\Delta t$ (ms)')
xlim([-15 15])
xticks([-15:5:15])
ylabel('Work ($MgL$)')
text(tmp,lblx,lbly,'F','units','normalized')

%legend('Ankle','Knee','Hip')
tmp2=tmp.YLim;
plot([0 0],tmp2,'color',[0 0 0 .5])
plot([delayVec(1) delayVec(end)],[0 0],'color',[0 0 0 .5])
yyaxis right
ylabel('(J)')
tmp.YColor=[0 0 0];
tmp.YLim=tmp2*Tnorm;
box off
%legend box off

% tmp=subplot(221);
% yyaxis left
% tmp.YColor=[0 0 0];
% plot(delayVec,ankle_angle_error,'k');hold on
% plot(delayVec,knee_angle_error,'k--');hold on
% plot(delayVec,hip_angle_error,'k:','linewidth',1.1);hold on
% text(tmp,lblx,lbly,'A','units','normalized')
% %xlabel('$\Delta t$ (ms)')
% ylabel('Joint angle (rad)')
% legend('Ankle','Knee','Hip')
% legend box off
% xlim([-15 15])
% xticks([-15:5:15])
% %title('Angle')
% tmp2=tmp.YLim;
% yyaxis right
% ylabel('(deg)')
% tmp.YColor=[0 0 0];
% tmp.YLim=tmp2*angle_norm;
% box off

tmp=subplot(322);
yyaxis left
tmp.YColor=[0 0 0];
plot(delayVec,LFT_T_prox_error/Tnorm,'k');hold on
plot(delayVec,LSK_T_prox_error/Tnorm,'k--');hold on
plot(delayVec,LTH_T_prox_error/Tnorm,'k:','linewidth',1.1);hold on
%legend('ankle','knee','hip','Location','north')
%legend box off
title('Walking data')
xlim([-15 15])
xticks([-15:5:15])
%xlabel('$\Delta t$ (ms)')
ylabel('Moment RMSE ($MgL$)') 
text(tmp,lblx,lbly,'D','units','normalized')
tmp2=tmp.YLim;
yyaxis right
ylabel('(N m)')
tmp.YColor=[0 0 0];
tmp.YLim=tmp2*Tnorm;
box off

tmp=subplot(324);
plot(delayVec,LeftAnklePower_error/Pnorm,'k');hold on
plot(delayVec,LeftKneePower_error/Pnorm,'k--');hold on
plot(delayVec,LeftHipPower_error/Pnorm,'k:','linewidth',1.1);hold on
%xlabel('$\Delta t$ (ms)')
%ylabel('($Mg^{1.5} L^{0.5}$)')
ylabel('Power RMSE ($Mg^{1.5} L^{0.5}$)')
text(tmp,lblx,lbly,'E','units','normalized')
xlim([-15 15])
xticks([-15:5:15])
tmp2=tmp.YLim;
yyaxis right
ylabel('(W)')
tmp.YColor=[0 0 0];
tmp.YLim=tmp2*Pnorm;
box off

T_sens=[LFT_T_prox_error(17) LSK_T_prox_error(17) LTH_T_prox_error(17)] %[Nm/ms]
T_sens_norm=T_sens/Tnorm/time_norm*1000 %[dimless]

P_sens=[LeftAnklePower_error(17) LeftKneePower_error(17) LeftHipPower_error(17)]
P_sens_norm=P_sens/Pnorm/time_norm*1000

T_sens_5=[LFT_T_prox_error(21) LSK_T_prox_error(21) LTH_T_prox_error(21)]
P_sens_5=[AnklePowerRMSE(21) LeftKneePower_error(21) LeftHipPower_error(21)]

AnklePowerRMSE_5_rel=AnklePowerRMSE_rel(21)
T_sens_5_rel=[LFT_T_prox_error_rel(21) LSK_T_prox_error_rel(21) LTH_T_prox_error_rel(21)]
HipPowerRMSE_5_rel=HipPowerRMSE_rel(21)

W_ankle_pos_0ms=W_ankle(16,1)/n_strides
W_ankle_pos_5ms=W_ankle(21,1)/n_strides
W_ankle_pos_per=W_ankle_pos_5ms/W_ankle_pos_0ms -1

W_hip_pos_0ms=W_hip(16,1)/n_strides
W_hip_pos_5ms=W_hip(21,1)/n_strides
W_hip_pos_per=W_hip_pos_5ms/W_hip_pos_0ms -1


%% GRF Fres and Mres
figure(23)
axlimit=[0 1.6 -400 400];
    % GRF
    tmp=subplot(3,4,1);plot(base,avgMod.grfl(:,2:3,idx_tmp)/Fnorm,'color',Lcolor); hold on
    plot(base,circshift(avgMod.grfr(:,2:3,idx_tmp)/Fnorm,0),'color',Rcolor); hold on
    yyaxis left
    tmp.YColor=[0 0 0];

    %subplot(3,3,iCond);plot(base,grf_tot(:,2:3,idx_tmp)/Fnorm,'color',totCol); hold on    
    ylim([-.5 1.5])
    ylabel('Force ($Mg$)')
    xlabel(TimeLabel)
    %text(tmp,lblx,lbly,'A','units','normalized')
    %legend('v=0.7 m/s','v=0.9 m/s','v=1.1 m/s')
    box off
    %title(titleName)
    tmp2=tmp.YLim;
    yyaxis right
    ylabel('(N)')
    tmp.YColor=[0 0 0];
    tmp.YLim=tmp2*Fnorm;
    box off

    
%     tmp=subplot(322);plot(base,avgMod.Msum_wrt_COM_lhs(:,3,idx_tmp)/Tnorm,'color',totCol); hold on
%     %plot(base,circshift(avgMod.Mfree_r(:,3,idx_tmp)/Tnorm,0),'color',rightCol); hold on
%     %subplot(3,3,iCond);plot(base,grf_tot(:,2:3,idx_tmp)/Fnorm,'color',totCol); hold on    
%     yyaxis left
%     tmp.YColor=[0 0 0];
%     text(tmp,lblx,lbly,'B','units','normalized')
%     ylim([-.1 .1])
%     ylabel('Ground reaction moment ($MgL$)')
%     xlabel(TimeLabel)
%     %legend('v=0.7 m/s','v=0.9 m/s','v=1.1 m/s')
%     tmp2=tmp.YLim;
%     yyaxis right
%     ylabel('(N m)')
%     tmp.YColor=[0 0 0];
%     tmp.YLim=tmp2*Tnorm;
%     box off
%     
    
    
    
     % FRES mag
     figure(69)
    tmp=subplot(2,4,7);plot(delayVec,mean_Fres,'color',Ncolor); hold on
    yyaxis left
    tmp.YColor=[0 0 0];
    %text(tmp,lblx,lbly,'G','units','normalized')
    ylabel({'Residual force';'RMSE ($Mg$)'})
    ylim([0 0.16])
    xlabel('Delay $\Delta t$ (ms)')
    xlim([-15 15])
    xticks([-15:5:15])
    tmp2=tmp.YLim;
    yyaxis right
    ylabel('(N)')
    tmp.YColor=[0 0 0];
    tmp.YLim=tmp2*Fnorm;
    box off
    
    % Mres
    tmp=subplot(2,4,8);plot(delayVec,mean_Mres,'color',Ncolor); hold on
    yyaxis left
    tmp.YColor=[0 0 0];
    ylim([0 0.05])
    ylabel({'Residual moment';'RMSE ($Mg$)'})
    xlabel('Delay $\Delta t$ (ms)')
    %text(tmp,lblx,lbly,'H','units','normalized')
    xlim([-15 15])
    xticks([-15:5:15])
    tmp2=tmp.YLim;
    yyaxis right
    ylabel('(N m)')
    tmp.YColor=[0 0 0];
    tmp.YLim=tmp2*Tnorm;
    box off   
    
    Fres_sens=diff(mean_Fres([25,31]))*Fnorm/6
    Fres_sens_norm=Fres_sens/Fnorm/time_norm*1000
    Mres_sens=diff(mean_Mres([25,31]))*Tnorm/6
    Mres_sens_norm=Mres_sens/Tnorm/time_norm*1000
    
    Fres_5=[mean_Fres(21)*Fnorm mean_Fres(21)/mean_Fres(16) mean_Fres(21)]
%% % NET WORK, RATE OF MECH ENERGY CHANGE and RESIDUAL
work_norm=n_strides_all;
sumJointWork_with_feet=sumJointWork_with_feet./work_norm;
sumJointWork_no_feet=sumJointWork_no_feet./work_norm;
W_com_PLUS_per=W_com_PLUS_per./work_norm;
mod.delta_Emech_model=mod.delta_Emech_model./work_norm;
Wsoft_with_feet=Wsoft_with_feet./work_norm;
Wres_net=Wres_net./work_norm;

% something like this to fix double axis:
% yyaxis left
% tmp.YColor=[0 0 0];
% ylabel('[norm]')
% yyaxis right
% tmp.YColor=[0 0 0];
% tmp.YLim=[-0.5 3]*force_norm;
% ylabel('[N]')
bar_time=0.2; % [s]
bar_length=bar_time./stride_time*100;

base=linspace(0,100,101);
airpumpCond=med;
Lcolor=[0.8 0 0];
Rcolor=[0 .49 .38];
Ncolor=[0 0 0];
alphaMod=linspace(0.3,1,length(idg_pref));alphaMod=fliplr(alphaMod);
TimeLabel=['Time [% gait cycle]'];
speed_labels=[{'0.7 m/s'}, {'0.9 m/s'}, {'1.1 m/s'}];

for itmp=1:length(stride_time)
    base_stride(itmp,:)=linspace(0,stride_time(itmp),101);
end



%% work fig

w_pos_air=46.81; % [J/stride]
w_neg_air=-7.65;
w_net_air=w_pos_air+w_neg_air;

fig=figure;
fig.Name=['all work terms'];
xLimit_work=[0.7 3.3];
W_ankle_tot=[sum(W_ankle(:,[1,3]),2) sum(W_ankle(:,[2,4]),2)];
W_knee_tot=[sum(W_knee(:,[1,3]),2) sum(W_knee(:,[2,4]),2)];
W_hip_tot=[sum(W_hip(:,[1,3]),2) sum(W_hip(:,[2,4]),2)];
W_fdist_tot=[sum(W_fdist(:,[1,3]),2) sum(W_fdist(:,[2,4]),2)];
thinCol=[0.8*Ncolor 0.5];
for iCond=1:3
    switch iCond
        case 1
            titleName=['normal preferred'];
            condition=idg_pref;
        case 2
            titleName=['normal constant sf'];
            condition=idg_csf;
        case 3
            titleName=['airpump preferred'];
            condition=airpumpCond;
    end
    % joint work
    subplot(3,3,iCond)
    for iSign=1:2
        % distal foot
        %plot(W_fdist(condition,iSign),'color',Lcolor); hold on
        %plot(W_fdist(condition,iSign),'v','color',Lcolor); hold on
        %plot(W_fdist(condition,iSign+2),'color',Rcolor); hold on
        %plot(W_fdist(condition,iSign+2),'v','color',Rcolor); hold on
        
        % ankle
        plot(belt_speeds(condition),W_ankle_tot(condition,iSign),'v','color',Ncolor,'MarkerSize',5); hold on
        plot(belt_speeds(condition),W_knee_tot(condition,iSign),'x','color',Ncolor,'MarkerSize',5); hold on
        plot(belt_speeds(condition),W_hip_tot(condition,iSign),'^','color',Ncolor,'MarkerSize',5); hold on
        
        plot(belt_speeds(condition),W_ankle_tot(condition,iSign),'color',thinCol); hold on
        %plot(W_ankle(condition,iSign+2),'color',Rcolor); hold on
        %plot(W_ankle(condition,iSign+2),'o','color',Rcolor); hold on
        
        % knee
        plot(belt_speeds(condition),W_knee_tot(condition,iSign),'color',thinCol); hold on
        %plot(W_knee(condition,iSign+2),'color',Rcolor); hold on
        %plot(W_knee(condition,iSign+2),'*','color',Rcolor); hold on
        
        % hip
        plot(belt_speeds(condition),W_hip_tot(condition,iSign),'color',thinCol); hold on
        %plot(W_hip(condition,iSign+2),'color',Rcolor); hold on
        %plot(W_hip(condition,iSign+2),'^','color',Rcolor); hold on
    end
    plot([0 belt_speeds(condition(end))+1],[0 0],'--','color',[0.8*Ncolor 0.7])
    box off
    ylim([-130 150])
    %xticklabels(speed_labels)
    title(titleName)
    switch iCond
        case 1
            ylabel('Work / stride [J]')
            xlim([0.5 2.4])
            xticks([.7 1.2 1.7 2.2])
        case 2
            legend('ankle','knee','hip','location','northwest');
            legend box off
            xlim([0.6 1.7])
            xticks([.7 1 1.3 1.6])
        case 3
            xlim([0.65 1.15])
            xticks([.7 0.9 1.1])
    end
    
    subplot(3,3,iCond+6)
    for iSign=1
        % ankle
        plot(belt_speeds(condition),sum(W_comPlusPer(condition,:),2),'+','color',Ncolor,'MarkerSize',5); hold on
        plot(belt_speeds(condition),sum(W_sumJoint_no_feet(condition,:),2),'>','color',Ncolor,'MarkerSize',5); hold on
        plot(belt_speeds(condition),sum(W_fdist_tot(condition,:),2),'s','color',Ncolor,'MarkerSize',5); hold on
        plot(belt_speeds(condition),sum(W_soft(condition,:),2),'o','color',Ncolor,'MarkerSize',5); hold on
        plot(belt_speeds(condition),coll_work_l(condition)+coll_work_r(condition),'*','color',Ncolor,'MarkerSize',5); hold on
        
        plot(belt_speeds(condition),coll_work_l(condition)+coll_work_r(condition),'color',thinCol); hold on
        
        plot(belt_speeds(condition),sum(W_fdist_tot(condition,:),2),'color',thinCol); hold on
        plot(belt_speeds(condition),sum(W_comPlusPer(condition,:),2),'color',thinCol); hold on
        
        plot(belt_speeds(condition),sum(W_sumJoint_no_feet(condition,:),2),'color',thinCol); hold on
        
        % hip
        plot(belt_speeds(condition),sum(W_soft(condition,:),2),'color',thinCol); hold on
        % airpump work
        if iCond==3
            plot([0.65 1.15],-2*[w_net_air w_net_air],'k--')
        end
    end
    plot([0.5/iCond belt_speeds(condition(end))+.3/iCond],[0 0],'--','color',[0.8*Ncolor 0.7])
    box off
    xlim([0.5 belt_speeds(condition(end))+.2])
    xlabel('treadmill speed')
    ylim([-130 150])
    
    switch iCond
        case 1
            ylabel('Net work / stride [J]')
            xlim([0.5 2.4])
            xticks([.7 1.2 1.7 2.2])
        case 2
            legend('COM+per','sum joints','soft tissue','distal foot','collision','location','north');
            legend box off
            xlim([0.6 1.7])
            xticks([.7 1 1.3 1.6])
        case 3
            xlim([0.65 1.15])
            xticks([.7 .9 1.1])
    end
    
    subplot(3,3,iCond+3)
    for iSign=1:2
        % ankle
        plot(belt_speeds(condition),W_comPlusPer(condition,iSign),'+','color',Ncolor,'MarkerSize',5); hold on
        plot(belt_speeds(condition),W_sumJoint_no_feet(condition,iSign),'>','color',Ncolor,'MarkerSize',5); hold on
        plot(belt_speeds(condition),W_fdist_tot(condition,iSign),'s','color',Ncolor,'MarkerSize',5); hold on
        plot(belt_speeds(condition),W_soft(condition,iSign),'o','color',Ncolor,'MarkerSize',5); hold on
        
        plot(belt_speeds(condition),W_fdist_tot(condition,iSign),'color',thinCol); hold on
        plot(belt_speeds(condition),W_comPlusPer(condition,iSign),'color',thinCol); hold on
        
        plot(belt_speeds(condition),W_sumJoint_no_feet(condition,iSign),'color',thinCol); hold on
        
        % hip
        plot(belt_speeds(condition),W_soft(condition,iSign),'color',thinCol); hold on
        
        if iCond==3
            switch iSign
                case 1
                    plot([0.65 1.15 ],-2*[w_net_air w_net_air],'k--')
                case 2
                    %plot([0.65 1.15 ],[w_neg_air w_neg_air],'k--')
            end
            
        end
        
    end
    plot([0.5/iCond belt_speeds(condition(end))+.3/iCond],[0 0],'--','color',[0.8*Ncolor 0.7])
    box off
    xlim([0.5 belt_speeds(condition(end))+.2])
    %xlabel('treadmill speed')
    ylim([-130 150])
    
    switch iCond
        case 1
            ylabel('Work / stride [J]')
            xlim([0.5 2.4])
            xticks([.7 1.2 1.7 2.2])
        case 2
            legend('COM+per','sum joints','soft tissue','distal foot','location','north');
            legend box off
            xlim([0.6 1.7])
            xticks([.7 1 1.3 1.6])
        case 3
            xlim([0.65 1.15])
            xticks([.7 .9 1.1])
    end
    % iCond
end % iSide
% aggregate work


%% individual leg GRF, Pcom, Pjoint, Pdist
fig=figure;
fig.Name=['individual leg GRF, Pcom, Pjoint, Pdist'];
yAxLim=[-450 450];
for iCond=1:3
    switch iCond
        case 1
            titleName=['normal preferred'];
            condition=idg_pref;
        case 2
            titleName=['normal constant sf'];
            condition=idg_csf;
        case 3
            titleName=['airpump preferred'];
            condition=airpumpCond;
    end
    for iSpeed=1:length(condition)
        leftCol=[Lcolor alphaMod(iSpeed)];
        rightCol=[Rcolor alphaMod(iSpeed)];
        totCol=[Ncolor alphaMod(iSpeed)];
        idx_tmp=condition(iSpeed);
        
        % GRF
        subplot(4,3,iCond);
        plot(base,avgMod.grfl(:,2:3,idx_tmp)/Fnorm,'color',leftCol); hold on
        plot(base,circshift(avgMod.grfr(:,2:3,idx_tmp)/Fnorm,50),'color',rightCol);
        ylim([-0.5 1.5])
        
        box off
        if iCond==1
            ylabel('GRF [BW]')
        end
        title(titleName)
        %legend('v=0.7 m/s','v=0.9 m/s','v=1.1 m/s')
        
        % Pcom
        subplot(4,3,iCond+3);
        plot(base,avgMod.p_coml(:,idx_tmp),'color',leftCol); hold on
        plot(base,circshift(avgMod.p_comr(:,idx_tmp),50),'color',rightCol);
        ylim([-0.4 1.5])
        if iCond==1
            ylabel('Pcom [W]')
        end
        %axis(axlimit)
        box off
        ylim([-500 500])
        % mod.Pper
        subplot(4,3,iCond+6);plot(base,avgMod.Pper(:,idx_tmp),'color',totCol); hold on
        if iCond==1
            ylabel('mod.Pper [W]')
        end
        box off
        ylim([-500 500])
        % Psoft
        subplot(4,3,iCond+9);plot(base,avgMod.P_com_PLUS_per(:,idx_tmp),'color',totCol); hold on
        if iCond==1
            ylabel('Pcom+per [W]')
        end
        xlabel(TimeLabel)
        %axis(axlimit)
        box off
    end
    ylim([-500 500])
end


%% abstract fig: soft tissue and body powers
if 0
    alphaMod=linspace(0.3,1,3);alphaMod=fliplr(alphaMod);
    fig=figure;
    fig.Name=['Abstract: total, distal foot and soft tissue power'];
    
    yAxLim=[-420 290];
    for iCond=[1,2]
        switch iCond
            case 1
                titleName=['normal'];
                condition=idg_pref;
            case 3
                titleName=['normal constant sf'];
                condition=idg_csf;
            case 2
                titleName=['artificial soft tissue'];
                condition=airpumpCond;
        end
        for iSpeed=1:3%:length(condition)
            leftCol=[Lcolor alphaMod(iSpeed)];
            rightCol=[Rcolor alphaMod(iSpeed)];
            totCol=[Ncolor alphaMod(iSpeed)];
            idx_tmp=condition(iSpeed);
            
            % Pcom + per
            subplot(4,2,iCond);
            plot(base,P_com_plus_per_avg(:,idx_tmp),'color',totCol); hold on
            plot(base,zeros(size(base)),'--','color',[Ncolor 0.3],'linewidth',0.4); hold on
            box off
            if iCond==1
                ylabel('$P_{\mathrm{com+per}} \, [\mathrm{W}]$','interpreter','latex','fontsize',10)
            end
            title(titleName)
            %legend('v=0.7 m/s','v=0.9 m/s','v=1.1 m/s')
            ylim(yAxLim)
            xticks([0 50 100])
            xticklabels([])
            
            % Pjoint
            subplot(4,2,iCond+2);plot(base,mod.sumJointPower_no_feet_avg(:,idx_tmp),'color',totCol); hold on
            plot(base,zeros(size(base)),'--','color',[Ncolor 0.3],'linewidth',0.4); hold on
            if iCond==1
                %ylabel('Sum Pjoint [W]')
                ylabel('$\Sigma P_{\mathrm{joint}} \, [\mathrm{W}]$','interpreter','latex','fontsize',10)
                
            end
            xticks([0 50 100])
            xticklabels([])
            %axis(axlimit)
            box off
            ylim(yAxLim)
            % Pfeet
            subplot(4,2,iCond+4);plot(base,Pground_lfoot_avg(:,idx_tmp)+circshift(Pground_rfoot_avg(:,idx_tmp),50),'color',totCol); hold on
            plot(base,zeros(size(base)),'--','color',[Ncolor 0.3],'linewidth',0.4); hold on
            if iCond==1
                %ylabel('Pdist [W]')
                ylabel('$P_{\mathrm{dist}} \, [\mathrm{W}]$','interpreter','latex','fontsize',10)
                
            end
            box off
            ylim(yAxLim)
            xticks([0 50 100])
            xticklabels([])
            
            % Psoft
            % plot the data
            subplot(4,2,iCond+6);plot(base,mod.Psoft_no_feet_avg(:,idx_tmp),'color',totCol); hold on
            % plot zero line
            plot(base,zeros(size(base)),'--','color',[Ncolor 0.3],'linewidth',0.4); hold on
            % plot patch to fill integral
            x_patch=[base fliplr(base)];
            y_patch=[zeros(size(base)) fliplr(mod.Psoft_no_feet_avg(:,idx_tmp)')];
            tmp_p=patch(x_patch,y_patch,Ncolor);
            tmp_p.EdgeAlpha=0;
            tmp_p.FaceAlpha=0.04;
            xlabel([])
            
            if iCond==1
                %ylabel('Psoft [W]')
                ylabel('$P_{\mathrm{soft}} \, [\mathrm{W}]$','interpreter','latex','fontsize',10)
                
            end
            xlabel('Time [\% gait cycle]','interpreter','latex','fontsize',10)
            %axis(axlimit)
            box off
            plot(7+[0 bar_length(idx_tmp)],-85*[4-iSpeed 4-iSpeed]+320,'color',totCol)
            text(24.5,-85*[4-iSpeed]+335,speed_labels(iSpeed),'fontsize',6,'interpreter','latex')
            if iSpeed==3
                text(0+bar_length(idx_tmp)/2,305,'0.2 s','fontsize',6,'interpreter','latex')
            end
        end
        ylim(yAxLim)
    end
    
    ha=get(fig,'children');
    for i_ha=1:length(ha)
        pos=ha(i_ha).Position;
        if i_ha==1 || i_ha==5
            %continue
        end
        pos(2)=pos(2)*.85;
        set(ha(i_ha),'position',pos)
        set(ha(i_ha),'TickLabelInterpreter','latex')
    end
end