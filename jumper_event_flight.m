function [ VALUE,ISTERMINAL,direction ] = jumper_event_flight( t,state,parms )
%[ VALUE,ISTERMINAL,DIRECTION ] = jumper_event( t,state,parms )
%   this event function detects when the jumper reaches the apex of the
%   flight phase (more or less) by checking when the upward velocity of the
%   base is zero
nseg=length(parms.segparms.L);

ISTERMINAL = 1; % stop integration
direction = -1; % terminate when velocity goes through zero coming from positive value

base=state(2*nseg+1:2*nseg+2);
based=state(2*nseg+3:2*nseg+4); % velocity of base
VALUE=based(2); % vertical velocity of base
%VALUE=base(2); % vertical position of base
end