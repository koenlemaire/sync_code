function [grfl, grfr, grml, grmr] = detrend_horizontal_grf (t,hsl,tol,hsr,tor,grfl,grfr,grml,grmr)
% makes the average horizontal grf_tot zero, while also making the grf hard
% zero during swing phase. 
idx_swing_left=[];
event_margin=3;
for i_step=1:length(tol)-1
    idx_swing_left=[idx_swing_left (tol(i_step)+event_margin):(hsl(i_step+1)-event_margin)];
end
idx_swing_right=[];
for i_step=1:length(tor)-1
    idx_swing_right=[idx_swing_right (tor(i_step)+event_margin):(hsr(i_step)-event_margin)];
end
% make forces and moments 0 during swing
idx_swing_left(idx_swing_left>length(grfl))=[];
idx_swing_right(idx_swing_right>length(grfr))=[];


grfl(idx_swing_left,:)=0;
grfr(idx_swing_right,:)=0;
grml(idx_swing_left,:)=0;
grmr(idx_swing_right,:)=0;

% compute average grf in horizontal plane:
grf_tot=grfl+grfr;
gr_impulse=trapz(grf_tot);
% stance phases:
idx_stance_left=setdiff(1:length(t),idx_swing_left);
idx_stance_right=setdiff(1:length(t),idx_swing_right);
% correct only stance phase:
grfl(idx_stance_left,1:2)=grfl(idx_stance_left,1:2)-...
    1/(length(idx_stance_left)-1)*gr_impulse(1:2)/2;
grfr(idx_stance_right,1:2)=grfr(idx_stance_right,1:2)-...
    1/(length(idx_stance_right)-1)*gr_impulse(1:2)/2;
end
