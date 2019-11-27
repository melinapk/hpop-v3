function  T = ECI2LVLH(Y0)


%-----------------------------------------------------------------------
%
% Converts Eath Centered Inertial (ECI) Coordinates to Local Coordinates
% (LVLH).
% INPUTS:
%
%   Y0(1:3) = [3 x 1] {column vector, numeric} position in ECI frame.
%
%   Y0(4:6) = [3 x 1] {column vector, numeric} velocity in the ECI frame.
%
% OUTPUTS:
%
%   T = [3 x 3] {array, numeric} coordinate transformation matrix for the
%       transformation from the Earth Centered Inertial frame to the Local-
%       Vertical/Local-Horizontal frame.
%   EULER ANGLES 
%    phi ;
%    theta ;
%    psi ; 
%   The aforementioned angles are going to be used for STK verification of
%   results. 
% 
%
% REFERENCE: Markley, F. Landis, and John L. Crassidis. Fundamentals of 
%       spacecraft attitude determination and control. Vol. 33. New York: 
%       Springer, 2014. Pg 36. Eq 2.78-79.
%------------------------------------------------------------------------

upveci = Y0/sqrt(sum(Y0.*Y0));

% Carculation of normal to orbit vectorn 
% Angular momentum

for i = 1:1:3
    r(i) = Y0(i);
    v(i) = Y0(i+3);
end


%CALULATION: LVLH X axis from position vector 

ver = - r/norm(r);

%CALCULATION: LVLH Y axis from velocity vector

t = v/norm(v);

%CALCULATION: LVLH Z axis from normal to orbit vector 

h = cross(t,ver);

%CTM from ECI to LVLH 

T = [t' , h' , ver'];

%-----------------------Phase A-------------------------------------
%       CALCULATE: Euler angles in degrees 
%-------------------------------------------------------------------


%q1 = rotm2quat(T)
seq = 'ZYX';
eul1 = rotm2eul(T, seq);
eul2 = rad2deg(eul1)

end