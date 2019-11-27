function G = input_matrix(orbit)

%% obtain orbital elements / position / angular momentum from orbit vector 

sm = orbit(1);
ecc = orbit(2);
RA = orbit(3);
inc = orbit(4);
ArgPer = orbit(5);
TrueA = orbit(6);
T = orbit(7);
b = orbit(8);
p = orbit(9);
r = orbit(10);
h = orbit(11);

%% calculate the input matrix of the first control system. - Orbit Control

G = [ (2*sm^2*ecc*sin(TrueA))/h , 2*sm^2*p/(h*r) , 0 ; p*sin(TrueA)/h , ((p+r)*cos(TrueA)+r*ecc)/h , 0 ; 0 , 0 , (r*cos(ArgPer+TrueA))/h ; 0 , 0 , r*sin(ArgPer+TrueA)/(h*sin(inc)) ; -(p*cos(TrueA))/(h*ecc) , ((p+r)*sin(TrueA))/(h*ecc) , -r*sin(ArgPer+TrueA)*cos(inc)/(h*sin(inc)) ];  
    