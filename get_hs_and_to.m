function [hsl, tol, hsr, tor] = get_hs_and_to (grfl_z,grfr_z,treshold,give_warning)
% finds heelstrikes and toe offs by means of a tresholding method, returns
% indices of heel strikes and toe offs such that all data starts and ends
% with a left heel strike

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

n_strides=length(hsl)-1;

% safety check on continuity of gait events:
if length(tor)==n_strides && ... % checks on nr of gait events
        length(hsr)==n_strides && ...
        length(tol)==n_strides % now we're good
    events_in_sequence=[hsl(1:n_strides)'; tor'; hsr'; tol'];
    events_in_sequence=events_in_sequence(:);
    if min(diff(events_in_sequence))<0 % some event is out of order
        isbad=1;
    else
        isbad=0;
    end
else
    isbad=1;
end

if max(diff(hsl))>1.5*mean(diff(hsl)) || ...
        max(diff(tol))>1.5*mean(diff(tol)) || ...
        max(diff(hsr))>1.5*mean(diff(hsr)) || ...
        max(diff(tor))>1.5*mean(diff(tor)) || ...
        min(diff(hsl))<0.5*mean(diff(hsl)) || ...
        min(diff(tol))<0.5*mean(diff(tol)) || ...
        min(diff(hsr))<0.5*mean(diff(hsr)) || ...
        min(diff(tor))<0.5*mean(diff(tor))
    isbad=1;
end

if give_warning && isbad    
    warning('there is probably a missing stride in this file, check fig and deal with it via keyboard')
%     fig=figure;
%     fig.Name=['heelstrike and toe off check'];
%     subplot(211);plot(grfl_z,'k'); hold on; plot(tol,grfl_z(tol,:),'bo'); plot(hsl,grfl_z(hsl,:),'ro')
%     legend('force','toe off','heel strike')
%     title('left')
%     subplot(212);plot(grfr_z,'k'); hold on; plot(tor,grfr_z(tor,:),'bo'); plot(hsr,grfr_z(hsr,:),'ro')
%     legend('force','toe off','heel strike')
%     title('right')
%     keyboard   
end