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

h= cross(r , v);   
hs = dot(h, h);   
hmag = sqrt(hs);

% unit angular momentum vector 

for i = 1:1:3    
    hhat(i) = h(i) / hmag;    
end

%Calsulation of position & velocity vector 

rv = dot(r,r);   
rmag = sqrt(rv);
vv = dot(v,v);
vmag = sqrt(vv);

% unit position vector

for i = 1:1:3   
    rhat(i) = Y0(i) / rmag;    
end

vhat = cross(hhat, rhat);


% LVLH axis unit vectors
o3 = -Y0(1:3)/norm(Y0(1:3));
o2 = -h/norm(h);
o1 = cross(o2,o3);



% CTM from ECI to LVLH
T = [o1', o2', o3'] ;


%-----------------------Phase A-------------------------------------
%       CALCULATE: Euler angles in degrees 
%-------------------------------------------------------------------



%q1 = rotm2quat(T)
% seq = 'ZYX';
% eul1 = rotm2eul(T, seq)
% 
% eul2 = rad2deg(eul1)
% r2 = eul2rotm(eul1,seq)
% q2 = rotm2quat(r2)



% %------------------------ Phase B ------------------------------------
% %       CALCULATE: initial Principle Rotations
% %---------------------------------------------------------------------
% 
% phi_i = atan2(eull(2,3),eullT(3,3));
% theta_i = -asin(eull(1,3));
% psi_i = atan2(eullT(1,2),eull(1,1));
% eul_ECI_i =[phi_i;theta_i;psi_i];
% 
% 
% e1 = phi_i;
% e2 = theta_i;
% e3 = psi_i;
% 
% 
% %-----------------------Phase C-------------------------------------
% %       CALCULATE: Euler angles in degrees 
% %-------------------------------------------------------------------
% 
% e1d = e1 * 180 / pi;
% e2d = e2 * 180 / pi;
% e3d = e3 * 180 / pi;

end