
function [sm, ecc, RA, Inclination, ArgPer, TrueA, T, b, p, r, h] = Orbital_Elem(Y)

%{
 This function computes the classical orbital elements (coe)
 from vector Y

mu - gravitational parameter (km^3/s^2)
R - position vector in the geocentric equatorial frame (km)
V - velocity vector in the geocentric equatorial frame (km)
r, v - the magnitudes of R and V
vr - radial velocity component (km/s)
H - the angular momentum vector (km^2/s)
h - the magnitude of H (km^2/s)
incl - inclination of the orbit (rad)
N - the node line vector (km^2/s)
n - the magnitude of N
cp - cross product of N and R
RA - right ascension of the ascending node (rad)
E - eccentricity vector
e - eccentricity (magnitude of E)
eps - a small number below which the eccentricity is considered
to be zero
w - argument of perigee (rad)
TA - true anomaly (rad)
sm - semimajor axis (km)
%}

mu = 398600.4; 

% i=0;

% for i = 1:1:1440 
%     

for k = 1:1:3 
%     R(k) = Eph(i, k+1)/1000;
%     V(k) = Eph(i, k+4)/1000;
R(k) = Y(k)/1000;
V(k) = Y(k+3)/1000;

end

eps = 1.e-10;

r = norm(R);
v = norm(V);
vr = dot(R,V)/r; 


H = cross(R,V);
h = norm(H);

%Calculate Inclination
incl = acos(H(3)/h);

Inclination =  rad2deg(incl);

%Calculate Node Line vector
ab = [ 0 , 0, 1];
N = cross(ab,H);
n = norm(N);


%Calculate right ascension 

if n ~= 0
Rasc = acos(N(1)/n);
if N(2) < 0
Rasc = 2*pi - Rasc;
end
else
Rasc = 0;
end

Ra = rad2deg(Rasc);
RA = Ra;

%Calculate eccentricity 
E = 1/mu*((v^2 - mu/r)*R - r*vr*V);
e = norm(E);

ecc = e;
%Calculate Argument of Perigee
if n ~= 0
if e > eps
w = acos(dot(N,E)/n/e);

if E(3) < 0
w = 2*pi - w;
end
else
w = 0;
end
else
w = 0;
end

W = rad2deg(w);

ArgPer = W ;

%Calculate true anomaly 
if e > eps
TA = acos(dot(E,R)/e/r);
if vr < 0
TA = 2*pi - TA;
end

else
    
cp = cross(N,R);
if cp(3) >= 0
TA = acos(dot(N,R)/n/r);
else
TA = 2*pi - acos(dot(N,R)/n/r);
end
end

TrueAnom = rad2deg(TA);
TrueA = TrueAnom;

%(sm < 0 for a hyperbola):
sm = h^2/mu/(1 - e^2); 

%Calculate period 
T = 2*pi/sqrt(mu)*sm^1.5;

% Calculate semi-minor axis
b = sm*sqrt(1-ecc^2);

% Calculate semi latus rectum 
p = sm*(1-ecc^2);


end 
