function  [Vector_A, Vector_B, Vector_C, Qxx , Plvlh ] = ECI2LVLHIV(Eph, vector1, vector2 )

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
%   The aforementioned angles are going to be used for STK verification of
%   results. 
%------------------------------------------------------------------------



for J = 1:1:1588 

% Carculation of normal to orbit vectorn 
% Angular momentum

for i = 1:1:3
    r(i) = Eph(J, i+1);
    v(i) = Eph(J, i+4);
end

i = 0;

%CALULATION: LVLH Z axis from position vector 

p = - r/norm(r);
for i = 1:1:3 
P(J,i) = p(i);
pv = dot(r,r);   
pmag = sqrt(pv);


%CONSTRUCTION OF STK'S VECTOR
%Creates the file format that is going to be imported in STK.

Vector_A(J,i+1) = P(J,i);
Vector_A(J,1) = Eph(J,1);

end

%CALCULATION: Velocity unit vector and magnitude 

i = 0;

t = v/norm(v);
for i = 1:1:3
V(J,i) = t(i);
tv = dot(v,v);
tmag = sqrt(tv);
end

i = 0;

%CALCULATION: LVLH Y axis from normal to orbit vector

hA = cross(r,v);
h = hA/norm(hA);
for i = 1:1:3
H(J,i) = h(i);
hs = dot(h, h);   
hmag = sqrt(hs);

%CONSTRUCTION OF STK'S VECTOR

Vector_B(J,i+1) = H(J,i);
Vector_B(J,1) = Eph(J,1);

end


i = 0;

%CALCULATIO: LVLH X axis from cross product of h x ver

cp = cross(h,p);
for i = 1:1:3
    
CP(J,i) = cp(i);
cps = dot(cp,cp);
cpmag = sqrt(cps);

%CONSTRUCTION OF STK'S VECTOR

Vector_C(J,i+1) = CP(J,i);
Vector_C(J,1) = Eph(J,1);

end

%-----------------------Phase A---------------------------------------
%       CALCULATE: Euler angles in degrees 
%	- The first rotation around the x axis of the I frame with ö angle
%	  rotates the the y axis of I frame to the y axis of the intermediate 
%	  frame I*
%	- The second rotation around the axis z of the I* frame, rotates 
%	  the x axis of the I* frame to the x axis of the LVLH frame
%---------------------------------------------------------------------


Qxx = [ -p ; cp ; h ];



% phirad = mod (real (sign(-cp(3))*acos(-cp(2)/sqrt(cp(2)^2+cp(3)^2))) , 2*pi);
% thetarad = acos(cp(1)/cpmag);
% 
% phi = rad2deg(phirad);
% theta = rad2deg(thetarad);
% psi = 0;
% 
% %CALCULATION: Elementary Rotation Matrices 
% 
% % R3 = [ cos(thetarad) , sin(thetarad) , 0 ; -sin(thetarad) , cos(thetarad) , 0 ; 0 , 0 , 1 ];
% % R2 = [ cos(psi) , 0 , -sin(psi) ; 0 , 1 , 0 ; sin(psi) , 0 , cos(psi) ];
% % R1 = [ 1, 0 , 0 ; 0 , cos(phirad) , sin(phirad) ; 0 , -sin(phirad) , cos(phirad) ];
% % 
% % 
% % RM = R3*R1;
% 
% Tx = [ 1 , 0 , 0 ; 0, cos(phi) , sin(phi) ; 0 , -sin(phi) , cos(phi) ];
% Ty = [ cos(theta) , 0 , -sin(theta) ; 0 , 1 , 0 ; sin(theta) , 0 , cos(theta) ];
% Tz = [ cos(psi) , sin(psi) , 0 ; -sin(psi) , cos(psi) , 0 ; 0 , 0 , 1 ]
% 
% TM = Tz*Ty*Tx ; 


for a = 1:1:3
    
    peci(a) = vector1(J,a);
    veci(a) = vector2(J,a);
    
end

plvlh = Qxx*peci';
vlvlh = Qxx*veci';

for k = 1:1:3 
    Plvlh(J,1) = Eph(J,1);
    Plvlh(J,k+1) = plvlh(k);
    
    Vlvlh(J,1) = Eph(J,1);
    Vlvlh(J,k+1) = vlvlh(k);
end

% seq = 'ZYX';
% eul1 = rotm2eul( Qxx , seq);
% eul2 = rad2deg(eul1);

alpha = atan2d_0_360(Qxx(3,1),-Qxx(3,2));
beta = acosd(Qxx(3,3));
gamma = atan2d_0_360(Qxx(1,3), Qxx(2,3));

   
Angles321(J,1) = Eph(J,1);
Angles321(J,2) = alpha;
Angles321(J,3) = beta;
Angles321(J,4) = gamma;

end

end