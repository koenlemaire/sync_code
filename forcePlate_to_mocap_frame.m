function [grfl,grfr,grml,grmr]=forcePlate_to_mocap_frame(grfl,grfr,grml,grmr,r_left_FP_origin_in_mocap_frame)
% [grfl,grfr,grml,grmr]=calibrateForces(grf_raw,grm_raw)
% calibrateForces takes the output generated from 'load(fileName)' (ie
% samples) and transforms this into left and right groundreaction forces
% (grf) and moments (grm) acting on on the subject, in a single reference
% frame that is attached to the belt (see below for definition). 
%
% force plate native coordinate system:
% ! Walking direction (forward) is such that the rollers that are driven by
% the motor (v-string) are in front of you!!  
% ! The positive directions are defined wrt to forces acting ON THE BELT,
% not on the subject
%
% left origin: left back corner of left plate
% right origin: right back corner of right plate
% both plates:
% +X to left
% +Y forward
% +Z down
% distance between left and right origin: 0.978 m
% forces in data should already be acting ON BELT, not on subject
%
% desired coordinate system (the one transformed to in this function):
    % origin: left back corner FP - r_left_FP_origin_in_mocap_frame
    % +X to right
    % +Z up
    % +Y forward
% DATA:
% x-direction distance between left and right FP origin: 0.978 m
% x-direction distance between left FP and mocap origin: 0.013 m     
% y-direction distance between FP and mocap origin: 1.623 m
% z-direction distance between FP and mocap origin: 0.015 m **
% 15 mm according to bertec documentation, uncertain if they reflect
% actual mocap origin ....   
% Koen Lemaire okt 2019

% translate right force plate to left origin:
r_trans=repmat([0.978 0 0],length(grmr),1); % [m] difference in origins
grmr=grmr+cross(r_trans,grfr); % add moment due to translation of origin

% now translate both forces to mocap origin:
r_trans=repmat(r_left_FP_origin_in_mocap_frame,length(grmr),1); % [m] difference in origins
grmr=grmr+cross(r_trans,grfr); % add moment due to translation of origin
grml=grml+cross(r_trans,grfl); % add moment due to translation of origin
