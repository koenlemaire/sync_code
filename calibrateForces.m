function [grfl,grfr,grml,grmr]=calibrateForces(grf_raw,grm_raw)
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
% forces in data are acting ON BELT, not on subject
%
% desired coordinate system (the one transformed to in this function):
% origin: left back corner of left plate
% +X to right
% +Y forward
% +Z up

% Koen Lemaire okt 2019

% minus sign in calibration matrix to get forces acting on subject:
forceCalibrationMatrix = -diag([500 500 1000]);
momentCalibrationMatrix = -diag([800 250 400]);

grfl=grf_raw(:,1:3)*forceCalibrationMatrix; % [N]
grml=grm_raw(:,1:3)*momentCalibrationMatrix; % [Nm]
grfr=grf_raw(:,4:6)*forceCalibrationMatrix; % [N]
grmr=grm_raw(:,4:6)*momentCalibrationMatrix; % [Nm]

% transform force data (no need to change origin here...):
R=[-1 0 0;
    0 1 0;
    0 0 -1]; % transformation based on above description
grfl=grfl*R;
grfr=grfr*R;

% transform moment data:
% first apply the same rotation as to force data (I think ;))
grml=grml*R;
grmr=grmr*R;

