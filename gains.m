function K = gains (orbit, dt)


%% gains : calculate the specific gains of the orbital element feedback problem. 

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

%% semi-major axis gain 

Ka = h^2/(4*sm^4*(1+ecc)^2)*(1/dt);

%% eccentricity gain 

Ke = h^2/(4*p^2)*(1/dt);

%% inclination gain 

Ki = ((h+ecc*h*cos(ArgPer + asin(ecc*sin(ArgPer))))/(p*(-1+ecc^2*sin(ArgPer)^2)))^2*(1/dt);

%% right ascension node gain 

Kra = ((h*sin(inc)*(-1 + ecc*sin(ArgPer+asin(ecc*cos(ArgPer)))))/(p*(1-ecc^2*cos(ArgPer)^2)))^2*(1/dt);

%% argument of periapse  gain 

Karg = (ecc^2*h^2)/(4*p^2)*(1 - ecc^2/4)*(1/dt);

%% mean anomaly gain 

Kma = (sm^2*ecc^2*h^2)/(4*b^2*p^2)*(1 - ecc^2/4)*(1/dt);

%% set gain vector 

K = [ 10^(-4) ; 0 ; 400 ; 0 ; 0];

%K = [Ka, Ke, Ki, Kra, Karg];



