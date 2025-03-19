function [stated, Vnew, succes] = segdyn (state, segparms, K, V, A_cons, b_cons)

% General implementation of the segment dynamics for a 2-D chain model of the skeleton
% consisting of nseg segments connected in frictionless hinge joints. 
% Richard Casius & Knoek van Soest, Fac. FBW, Vrije Universiteit Amsterdam.
% 
% Simplified version: - requires exactly 3*nseg unknowns (one more for each constraint
%                       equation), i.e., least squares solutions are not supported
%                     - the complete system of equations is set up each time this
%                       function is invoked, i.e., very inefficient
% INPUT PARAMETERS
%
%    state:    vector containing 2*nseg+4 state variables: 
%                    [phi_1 .. phi_n, phid_1 .. phid_n, x, y, Vx, Vy], where 
%              phi_i is ith segment angle, phid_i the ith segment angular velocity, 
%              (x, y) the base position, and (Vx, Vy) the base velocity
%
%    segparms: vector of length 4*nseg: 
%                    [l_1 .. l_n, d_1 .. d_n, m_1 .. m_n, j_1..j_n], where
%              l_i is the length of the ith segment (m), d_i is the distance from the
%              proximal end of the ith segment to its center of gravity (m), m_i is the
%              mass of the ith segment (kg), and j_i the moment of inertia of the ith
%              segment with respect to its center of gravity (Nm/kg); Note: segparms may
%              also be a struct containing the vectors L, d, m and j, where 
%              L=[l_1 .. l_n], d=[d_1 .. d_n], m=[m_1 .. m_n], and j=[j_1 .. j_n]
%
%    K:        vector containg 7*nseg+5 truth values; if K(i)=1, V(i) should contain 
%              an appropriate value; if K(i)=0, V(i) is considered unknown
%
%    V:        vector containg 7*nseg+5 values for the follwing variables:
%                  Fxr   = V(1    : n+1);   nseg+1 horizontal joint reaction forces
%                  Fyr   = V(n+2  : n2+2);  nseg+1 vertical joint reaction forces
%                  M     = V(n2+3 : n3+3);  nseg+1 net joint moments
%                  Fxext = V(n3+4 : n4+3);  nseg horizontal external forces
%                  Fyext = V(n4+4 : n5+3);  nseg vertical externalforces
%                  Mext  = V(n5+4 : n6+3);  nseg external moments
%                  phidd = V(n6+4 : n7+3);  nseg segment angular accelerations
%                  Ax    = V(n7+4);         horizontal base acceleration
%                  Ay    = V(n7+5);         vertical base acceleration
%              exactly 4*nseg+5 of these variables should contain an appropriate value
%              (these 'known' variables must be identified through the array K)
%
%    A_cons:   left-hand side of constraint equations (optional)
%    b_cons:   right-hand side of constraint equations (optional)
%  
% OUTPUT PARAMETERS
%
%    stated:   first order derivative of state
%    Vnew:     copy of V, with all unknown variables replaced with the value resulting
%              from solving the equations of motion
%    succes:   1 if the equations of motion were solved successfully, 0 otherwise
 
first  = 1;     % may not be changed; just to indicate which parts need only be
                % exectuted once if all relevant information is stored
slow   = 0;     % may be changed; just to make Matlab-specific code more human readable
succes = 0;     % will be set to 1 if the equations of motion are solved successfully
stated = state; % get the right dimensions for ‘stated’ before the input vectors
                % are turned into columns vectors   
Vnew   = V;     % idem for Vnew

% make column vectors
state    = state(:);
K        = K(:);
V        = V(:);
if (first)
    if isa(segparms, 'double')
        segparms = segparms(:);
    elseif isa(segparms,'struct')
        if isfield(segparms,'L')
            l=segparms.L;
        elseif isfield(segparms,'l')
            l=segparms.l;
        else
            disp ('segparms does not contain field L');
            return;
        end
        if isfield(segparms,'d')
            d=segparms.d;
        elseif isfield(segparms,'D')
            d=segparms.D;
        else
            disp ('segparms does not contain field d');
            return;
        end
        if isfield(segparms,'m')
            m=segparms.m;
        elseif isfield(segparms,'M')
            m=segparms.M;
        else
            disp ('segparms does not contain field m');
            return;
        end
        if isfield(segparms,'j')
            J=segparms.j;
        elseif isfield(segparms,'J')
            J=segparms.J;
        else
            disp ('segparms does not contain field j');
            return;
        end
        segparms = [l(:) ; d(:) ; m(:) ; J(:)];  
    else
        disp ('segparms is of unsupported type; should be either a vector or a struct');        
    end
end 
 
% very basic input checks
if (first)  
   lk = length(K);
   if (length(V) ~= lk)
      disp ('Size of K and V differ.');
      return;
   end
   n  = (lk-5) / 7;
   ls = length (segparms); 
   if (ls ~= 4*n)
      err = sprintf('Length of K and V (%d) disagrees with length segparms (%d).',lk,ls); 
      disp (err);
      disp ('K and V should should contain 7*nseg+5 elements, and segparms 4*nseg.');
      return;
   end   
   ls = length(state); 
   if (ls ~= 2*n+4)
      err = sprintf ('Length of K and V (%d) disagrees with length state (%d).', lk,ls); 
      disp (err);
      disp ('K and V should should contain 7*nseg+5 elements, and state 2*nseg+4.');
      return;
   end   
   fK = find(K);
   fV = find(isfinite(V));
   if ((length(fK) ~= length(fV)) | (fK ~= fV))
      % K and V disagree; report first illegal element
      for i=1:7*n+5
         if (K(i) == 1)
            if (~isfinite(V(i)))
               err = sprintf ('Element %d in K and V does not agree.', i);
               disp(err);
               return;
            end
         elseif (K(i) ~= 0)
            err = sprintf ('Element %d in K should be either 0 or 1.', i);
            disp (err);  
            return;
         end
      end
   end    
   lu = length(find(~K));
   if (nargin > 4)
       lc = length(A_cons);
   else
       lc = 0;
   end
   if (lu > 3*n + lc)
      err=sprintf('%d unknowns found: at most %d elements in K may be zero!',lu,3*n+lc); 
      disp(err);
      return;
   end    
end % input checks succeeded
 
if (first)
   % auxiliary variables
   n2 = n+n;
   n3 = n+n2; 
   n4 = n+n3;
   n5 = n+n4;
   n6 = n+n5;
   n7 = n+n6;
   
   % extract segment parameters
   l = segparms(1    : n);
   d = segparms(n+1  : n2);
   m = segparms(n2+1 : n3);
   J = segparms(n3+1 : n4);
end

% auxiliary variables that depend on the current state variables
phi  = state(1:n);
phid = state(n+1:n2);
if (slow)
   p = l-d;
   s = sin(phi);
   c = cos(phi);
   for i=1:n
      phidsqr(i) = phid(i)^2;
      ds(i)      = d(i)*s(i);
      dc(i)      = d(i)*c(i);
      ps(i)      = p(i)*s(i);
      pc(i)      = p(i)*c(i);
      ls(i)      = l(i)*s(i);
      lc(i)      = l(i)*c(i);
   end
else 
   p       = l-d;
   phidsqr = phid.^2;
   diagd   = diag(d);
   diagp   = diag(p);
   diagl   = diag(l);
   diags   = diag(sin(phi));
   diagc   = diag(cos(phi));
   ds      = diag(diagd*diags);
   dc      = diag(diagd*diagc);
   ps      = diag(diagp*diags);
   pc      = diag(diagp*diagc);
   ls      = diag(diagl*diags);
   lc      = diag(diagl*diagc);
end
 
% below the main algorithm is listed; it closely follows the discourse from section 3.3 

% build static blocks
if (first)

   % Blocks 1,2,3,4,6 
   B1 = zeros(n2, n2+2);
   B2 = zeros(n2, n+1);
   B3 = diag(ones(n2,1));
   B4 = zeros(n2, n);
   B6 = zeros(n2, 2);
   j  = 1;
   for i=1:n2   
      B1(i,j:j+1) = [1 -1];
      if (i==n)
         j = j+2;
      else
         j = j+1;
      end
   end
   B6(1:n,1)    = -m;
   B6(n+1:n2,2) = -m;

   % Blocks 8,9,10,11,12
   B8  = zeros(n, n+1);
   B9  = zeros(n, n2);
   B10 = diag(ones(n,1));
   B11 = diag(-J);
   B12 = zeros(n, 2);
   j   = 1;
   for i=1:n   
      B8(i,j:j+1) = [1 -1];
      j = j+1;
   end
   
   % initialize Blocks 5 and 7
   B5 = zeros(n2, n);
   B7 = zeros(n, n2+2);
end

% Block 5
for i=1:n
   for j=1:i-1
      B5(i,j)   =  m(i)*ls(j);
      B5(i+n,j) = -m(i)*lc(j);
   end   
   B5(i,i)   =  m(i)*ds(i);
   B5(i+n,i) = -m(i)*dc(i);
end

% Block 7
for i=1:n
    B7(i,i:i+1)       = [ ds(i)  ps(i)];
    B7(i,i+n+1:i+n+2) = [-dc(i) -pc(i)];
end   
Atmp = [[B1;B7] [B2;B8] [B3;B9] [B4;B10] [B5;B11] [B6;B12]];
% initial right-hand side
b = zeros(n3,1);
for i=1:n
   for j=1:i-1
      b(i)   = b(i)   - m(i)*lc(j)*phidsqr(j);
      b(i+n) = b(i+n) - m(i)*ls(j)*phidsqr(j);
   end  
   b(i)   = b(i)   - m(i)*dc(i) * phidsqr(i);
   b(i+n) = b(i+n) - m(i)*ds(i) * phidsqr(i);
end   

if nargin > 5 % only if A_cons and b_cons are specified
    Atmp = [Atmp; A_cons];
    b    = [b;    b_cons];
end

% move known parts to the right hand side 
nrows = length(b); % may be larger than 3*nseg if constraints are added at some time
A     = zeros (nrows,nrows);
j     = 1;
for i=1:n7+5
   if (K(i))
      b = b - V(i)*Atmp(:,i);
   else   
      A(:,j) = Atmp(:,i);
      j = j+1;
   end
end   

x=A\b; % now solve the square system

% return the solution
if (length(find(~isfinite(x))) == 0)
   if (slow)
      j = 1;
      for i=1:n7+5
         if (~K(i))
            Vnew(i) = x(j);
            j = j+1;
         end
      end
   else
      Vnew(find(~K)) = x;
   end
   stated(1:n)    = phid;            % phid
   stated(n+1:n2) = Vnew(n6+4:n7+3); % phidd
   stated(n2+1)   = state(n2+3);     % Vx
   stated(n2+2)   = state(n2+4);     % Vy
   stated(n2+3)   = Vnew(n7+4);      % Ax
   stated(n2+4)   = Vnew(n7+5);      % Ay
   succes = 1;
else
   disp ('Equations of motion could not be solved.');
end
