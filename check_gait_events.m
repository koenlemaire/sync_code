function [isbad] = check_gait_events (hsl, tol, hsr, tor, grfl_z, grfr_z, do_figure)
% this function checks for the order and consistency of gait events:
% heelstrike left (hsl) - toe off right (tor) - heel strike right - toe off left
% furthermore it does some basic checks on the distance between events and
% the total nr of events. This function assumes the output of
% get_gait_events
% 
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

if do_figure || isbad    
    warning('there is probably a missing stride in this file, check fig')
    fig=figure;
    fig.Name=['heelstrike and toe off check'];
    subplot(211);plot(grfl_z,'k'); hold on; plot(tol,grfl_z(tol,:),'bo'); plot(hsl,grfl_z(hsl,:),'ro')
    legend('force','toe off','heel strike')
    title('left')
    subplot(212);plot(grfr_z,'k'); hold on; plot(tor,grfr_z(tor,:),'bo'); plot(hsr,grfr_z(hsr,:),'ro')
    legend('force','toe off','heel strike')
    title('right') 
end