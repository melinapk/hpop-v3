function  Qxx = LVLHTransformation(Y)

% N = size(Eph,1);

% for J = 1:1:N 

% Carculation of normal to orbit vectorn 
% Angular momentum

for i = 1:1:3
    r(i) = Y(i);
    v(i) = Y(i+3);
end

i = 0;

%CALULATION: LVLH Z axis from position vector 

x = r/norm(r);
% for i = 1:1:3 
% P(N,i) = x(i);
% pv = dot(r,r);   
% pmag = sqrt(pv);
% end

%CALCULATION: Velocity unit vector and magnitude 

i = 0;

t = v/norm(v);
% for i = 1:1:3
% V(N,i) = t(i);
% tv = dot(v,v);
% tmag = sqrt(tv);
% end

i = 0;

%CALCULATION: LVLH Y axis from normal to orbit vector

hA = cross(r,v);
z = hA/norm(hA);
% for i = 1:1:3
% H(N,i) = z(i);
% hs = dot(z, z);   
% hmag = sqrt(hs);
% end

i = 0;

%CALCULATIO: LVLH X axis from cross product of h x ver

cp = cross(hA,r);
y = cp/norm(cp);
% for i = 1:1:3
%     
% CP(N,i) = cp(i);
% cps = dot(cp,cp);
% cpmag = sqrt(cps);
% 
% end

Qxx = [ x ; y ; z ];

% for a = 1:1:3
%     
%     peci(a) = vector1(J,a);
%     veci(a) = vector2(J,a);
%     
% end
% 
% plvlh = Qxx*peci';
% vlvlh = Qxx*veci';


% end
end



