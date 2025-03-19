function Pjo_trans = CalcTransJointPower(omega_s1,v_com_s1,r_com_s1,v_prox_s2,r_prox_s2,F_s1_ON_s2)
% Pjo_trans = CalcTransJointPower(omega_s1,v_com_s1,r_com_s1,v_prox_s2,r_prox_s2,F_s1_ON_s2)
% Calculates translational joint power in an articulation between two
% segments s1 and s2. The net joint force acts from s1 on s2 (F_s1_ON_s2)
% at the "joint centre", the location of which is assumed to coincide with
% the proximal end of s2 (r_prox_s2). The underlying model is that s1 and
% s2 are both rigid bodies, which are connected by a massless 6-DOF actuator
% representing the 'joint'. It turns out that for such a model the net
% power exerted by the actuator can be calculated as follows:
%
% Procedure:
% The translational joint power is the sum of the translational power
% acting from s1 on s2 and that acting from s2 on s1. This is 
% calculated as the dot product of the joint contact force with the
% difference between the velocity of the joint centre when that is taken to
% be rigidly connected to s2 and when that is taken to be rigidly connected
% to s1. Note that because of this, the power is not sensitive to the
% choice of reference frame. Also note that the appearance of the 'velocity
% of the joint centre when that is taken to be rigidly connected to s1' is
% only a velocity by happenstance; it turns out that if the reaction force
% is assumed to act at a location P on s1, then the mechanics conspire to
% result in this apparent velocity. Contact Koen Lemaire
% (koenlemaire@posteo.net) for detailed explanation. 
%
% !! All vectors MUST be expressed in the same (eg global) reference frame
%
% definition of inputs: 
% r_... [n x 3] position vector
% v_... [n x 3] velocity vector 
% omega_s1 [n x 3] angular velocity vector of s1
% F_s1_ON_s2 [n x 3] force from s1 on s2 acting at proximal end of s2
%
% ! Note that the calculation presented here only captures the net power
% exerted by the massless joint actuator. The power exerted by the actuator
% on s1 and s2 depends on the choice of reference frame. 

% Koen Lemaire 6/2020

% vector from s1 centre of mass (r_com_s1) to joint centre (r_prox_s2)
r=r_prox_s2-r_com_s1;

% velocity of joint centre on s1 
v_jc_s1=v_com_s1 + cross(omega_s1,r,2);

% velocity of joint centre on s2 
v_jc_s2=v_prox_s2;

% power of joint contact force:
Pjo_trans=sum((v_jc_s2-v_jc_s1).*F_s1_ON_s2,2);
