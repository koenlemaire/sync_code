function [p_coml, p_comr, v_com_belt, grfl_new, grfr_new, varargout] = calc_com_power(grfl,grfr,v_belt,t,do_xyForce_correction)
% [p_coml_belt, p_comr_belt, v_com_belt, a_com_lab v_com_lab, subj_mass, w_com_belt, dEkin_com_belt dv_com_lab dEkin_com_belt_uncorrected dv_com_lab_uncorrected] = calc_com_power(grfl,grfr,v_belt,t)
% 
% calc_com_power calculates left and right leg contributions to com power
% (i.e. the rate of change of kinetic and potential gravitational energy of
% the centre of mass, in a belt fixed reference frame) by inner product of
% left and right leg ground reaction force with centre of mass velocity
% (v_com). v_com is calculated based on the total ground reaction force and
% expressed in a belt fixed reference frame such that mean(v_com) = [0
% v_avg 0]. Subj_mass is derived from the mean vertical grf. Compared to
% previous work, this function implements an improved way of dealing with
% trends in the horizontal ground reaction force, such that:
% 1) com velocity is very near periodic
% 2) net com work is consistent with the mechanical energy change of the
% com.   
%
% USED REFERENCE FRAME [X Y Z]:
%  +X right (wrt walking direction)
%  +Y forward (walking direction)
%  +Z upward
%
% INPUTS: 
%       grfl/grfr   [nx3], [N], left/right ground reaction force
%       v_belt      [scalar], [m/s], average forward velocity
%       t           [nx1], [s], time base for grf
%
% OUTPUTS:
%       p_coml/p_comr   [nx1], [W], left/right com power in treadmill frame
%       v_com           [nx1], [m/s], com velocity in treadmill/lab frame
%       subj_mass       [scalar], [kg], subject mass 
%       w_com*          [scalar], [J], net com work in treadmill frame
%       dEkin_com*      [scalar], [J], change in com kinetic energy in
%                           treadmill frame
%
% ! This function forces close to periodic motion, and will violate the
% data if the original data is not close to periodic. The user is
% encouraged to input 1) a large nr of strides and 2) an integer nr of
% strides. 
% ! This function supports non-equally spaced time bases; all calculations
% of mean values are based on trapz.
%
% Note on previous work (function calculateComWorkRateForSteadyStateGait
% from Karl Zelik): therein the correction for horizontal drift was applied
% to the accelerations but not to the forces, leading to inconsistency
% between com work rate and com mechanical energy. That issue is resolved
% here.

% Koen Lemaire 2020, inspired by Karl Zelik (2009)

% Joint power is independent of reference frame because it is an internal
% power. DeltaEkin = sum(Wint + Wext) and the dependency of DeltaEkin on
% reference frame arises through Wext .... The same goes for peripheral
% work/power. Here Wext is defined as work done by force acting on the
% (nonrigid!) body (gravity and grf) and Wint is defined as work done
% between segments (ie joint work). 
%% Subject mass
T=t(end)-t(1); % total time
g=-9.81; % [m/s^2] gravitational acceleration 
grf_z=grfl(:,3)+grfr(:,3); % [N] upward ground reaction force
subj_mass=-trapz(t,grf_z)/T/g; % [kg] mass of the subject

%% uncorrected change in com velocity and kinetic energy (sanity check):
f_sum=[grfl(:,1:2)+grfr(:,1:2) grf_z+g*subj_mass]; % [N] sum of forces on subject
P_sum_uncorrected=trapz(t,f_sum); % [Ns] net impulse on subject, uncorrected
vfinal_com_belt_uncorrected=P_sum_uncorrected/subj_mass+[0 v_belt 0]; % [m/s] final com velocity, in belt frame 
dEkin_com_belt_uncorrected=.5*subj_mass*(vfinal_com_belt_uncorrected*vfinal_com_belt_uncorrected' -v_belt.^2); % [J] change in COM kinetic energy, in belt frame

%% Fix trends in horizontal forces
% we are going to distribute equally the total trendline over
% the left and right leg.
% remove bias
F_avg_left=trapz(t,grfl(:,1:2))/T; % [N]
F_avg_right=trapz(t,grfr(:,1:2))/T; % [N]
F_avg_xy=F_avg_left+F_avg_right; % [N]

if do_xyForce_correction
    grfl(:,1:2)=grfl(:,1:2)-.5*F_avg_xy; % [N]
    grfr(:,1:2)=grfr(:,1:2)-.5*F_avg_xy; % [N]
end

grfl_new=grfl;
grfr_new=grfr;

% remove trendline (leads to more numerical error in velocity)
%trendline=tmp(:,1:2)-detrend(tmp(:,1:2)); % [N] the linear trend in the xy plane
%grfl(:,1:2)=grfl(:,1:2)-.5*trendline; % [N]
%grfr(:,1:2)=grfr(:,1:2)-.5*trendline; % [N]
%% Compute (corrected) COM velocity
f_sum=[grfl_new(:,1:2)+grfr_new(:,1:2) grf_z+g*subj_mass]; % [N] sum of forces on subject
a_com_lab=f_sum/subj_mass; % mean zero in all three dimensions! 
v_com_lab=cumtrapz(t,a_com_lab); % [m/s] delta com velocity in world frame, periodic because acc is mean zero ... 
v_com_lab=v_com_lab-trapz(t,v_com_lab)/T; % [m/s] integration constant such that mean velocity is zero in all directions
v_com_belt=[v_com_lab(:,1) v_com_lab(:,2)+v_belt  v_com_lab(:,3)]; % [m/s] transform com velocity to treadmill frame

%% Compute COM power
% COM power (in treadmill frame), grf fully consistent with v_com .... 
p_coml=dot(grfl_new,v_com_belt,2); % [W] left 
p_comr=dot(grfr_new,v_com_belt,2); % [W] right

%% Additional sanity checks (optional outputs)
% net com work in treadmill frame:
w_com_belt=trapz(t,p_coml+p_comr); % [J]
% change in com kinetic energy in treadmill frame:
dEkin_com_belt=.5*subj_mass*(v_com_belt(end,:)*v_com_belt(end,:)' - v_com_belt(1,:)*v_com_belt(1,:)'); % [J]
% change in com velocity (check on periodicity):
dv_com_lab=v_com_lab(end,:)-v_com_lab(1,:); % [m/s]

% [p_coml_belt, p_comr_belt, v_com_belt, v_com_lab, subj_mass, w_com_belt, dEkin_com_belt dv_com_lab dEkin_com_belt_uncorrected dv_com_lab_uncorrected] = calc_com_power(grfl,grfr,v_belt,t)
optout=[{a_com_lab},{v_com_lab}, {subj_mass}, {w_com_belt}, {dEkin_com_belt}, ...
    {dv_com_lab}, {dEkin_com_belt_uncorrected}, {F_avg_xy}, {F_avg_left}, {F_avg_right}];
n_out=nargout;
if n_out > 5
    for i=1:n_out-5
        varargout{i}=optout{i};
    end
end