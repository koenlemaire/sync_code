function Pjo_6DOF = Calc_6DOFJointPower(omega_s1,omega_s2,v_com_s1,r_com_s1,v_prox_s2,r_prox_s2,F_s1_ON_s2,T_s1_ON_s2)
% Pjo_6DOF = Calc_6DOFJointPower(omega_s1,omega_s2,v_com_s1,r_com_s1,v_prox_s2,r_prox_s2,F_s1_ON_s2,T_s1_ON_s2)
%
% Calculates joint power in a 6-Degree-Of-Freedom articulation between two
% segments s1 and s2. The net joint force acts from s1 on s2 (F_s1_ON_s2)
% at the "joint centre", the location of which is assumed to coincide with
% the proximal end of s2 (r_prox_s2). The underlying model is that s1 and
% s2 are both rigid bodies, which are connected by a massless 6-DOF actuator
% representing the 'joint'. Note that the force/torque from s1 on s2 are
% really the force/torque from the actuator on s2, the torque acting on s1
% depends on the location of the distal end of s1, but it turns out that
% the calculation of net joint power is independent of this locaiton... For
% this model the net power exerted by the actuator can be calculated as
% follows.
%
% Procedure:
% The 6-DOF joint power is the sum of the power acting from the actuator on
% s2 and that acting from the actuator on s1. Taking the massless 
% actuator as the free body, it can be shown that this joint power equals
% the the dot product of the joint contact force with the difference
% between the velocity of the joint centre when that is taken to be rigidly
% connected to s2 and when that is taken to be rigidly connected 
% to s1, plus the dot product of the joint angular velocity with the joint
% torque. Note that the distinction between rotational and translational
% power in reality does not make sense; part of the 'translational' power
% is in fact due to a torque acting within the massless actuator, ie the
% 'pure' torque acting from the actuator on s1 is not equal to that acting
% from the actuator on s2 ... In fact, calculating any '3-DOF' power when
% the dynamics allow for a 6-DOF joint plainly does not make sense ...
% Contact Koen Lemaire (koenlemaire@posteo.net) for detailed explanation. 
%
% Note that the joint power is by definition (and in this calculation) not
% sensitive to the choice of reference frame. 
%
% !! All vectors MUST be expressed in the same (eg global) reference frame
%
% definition of inputs: 
% r_... [n x 3] position vector
% v_... [n x 3] velocity vector 
% omega_... [n x 3] angular velocity vector of s1/s2
% F/T_s1_ON_s2 [n x 3] force/torque from s1 (or actually from massless
% actuator) on s2 acting at proximal end of s2 
%
% ! Note that the calculation presented here only captures the net power
% exerted by the massless joint actuator. The power exerted by the actuator
% on s1/s2 depends on the choice of reference frame. 

% Koen Lemaire 6/2020

% vector from s1 centre of mass (r_com_s1) to joint centre (r_prox_s2)
r=r_prox_s2-r_com_s1;

% velocity of joint centre on s1 
v_jc_s1=v_com_s1 + cross(omega_s1,r,2);

% velocity of joint centre on s2 
v_jc_s2=v_prox_s2;

% "power of joint contact force":
Pjo_trans=sum((v_jc_s2-v_jc_s1).*F_s1_ON_s2,2);

% "rotational joint power":
Pjo_rot=sum((omega_s2-omega_s1).*T_s1_ON_s2,2);

Pjo_6DOF=Pjo_trans+Pjo_rot; % note that this partition is a bit misleading
