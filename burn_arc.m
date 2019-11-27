function t_burn = burn_arc(orbit)

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

mu = 398600.4; 

cosE = (ecc + cos(TrueA))/(1+ecc*cos(TrueA));

E = acos (cosE); 

t_burn = 2*sqrt(sm^3/mu)*(E + ecc*sin(E));