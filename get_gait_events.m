function [hsl, tol, hsr, tor] = get_gait_events (grfl_z,grfr_z,treshold,give_warning)
% finds heelstrikes and toe offs by means of a tresholding method, returns
% indices of heel strikes and toe offs such that all data starts and ends
% with a left heel strike. Call check_gait_events to check the results of
% this function!

mpd=30; % samples, minimum distance between events
sample_margin=5; % nr samples away from the detection treshold

% left foot heel strike
idx=find(grfl_z<treshold); % sample indices where foot is off the treadmill
[~,idx2]=findpeaks(diff(idx),'minpeakdistance',mpd); % samples in idx where there's a 'jump'
hsl=idx(idx2)-sample_margin; % samples indices where heel strike happens

% left foot toe off
idx=find(grfl_z>=treshold); % sample indices where foot is on the treadmill
[~,idx2]=findpeaks(diff(idx),'minpeakdistance',mpd); % samples in idx where there's a 'jump'
tol=idx(idx2)+sample_margin; % samples indices where heel strike happens

% right foot heel strike
idx=find(grfr_z<treshold); % sample indices where foot is off the treadmill
[~,idx2]=findpeaks(diff(idx),'minpeakdistance',mpd); % samples in idx where there's a 'jump'
hsr=idx(idx2)-sample_margin; % samples indices where heel strike happens

% right foot toe off
idx=find(grfr_z>=treshold); % sample indices where foot is off the treadmill
[~,idx2]=findpeaks(diff(idx),'minpeakdistance',mpd); % samples in idx where there's a 'jump'
tor=idx(idx2)+sample_margin; % samples indices where heel strike happens

% make everything start with first left heel strike
tol(tol>hsl(end)|tol<hsl(1))=[];
hsr(hsr>hsl(end)|hsr<hsl(1))=[];
tor(tor>hsl(end)|tor<hsl(1))=[];

