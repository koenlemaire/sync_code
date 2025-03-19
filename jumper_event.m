function [ VALUE,ISTERMINAL,direction ] = jumper_event( t,state,parms )
%[ VALUE,ISTERMINAL,DIRECTION ] = jumper_event( t,state,parms )
%   this event function detects when the jumper leaves the ground by
%   checking when the y_grf is zero
[statedot,output]=frwd_jumper(t,state,parms); % unfortunately we have to run this again ... 


ISTERMINAL = 1;
direction = -1;

VALUE=output(6); % vertical component of grf
end