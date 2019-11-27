function Kerr = error_function(orbit)

global master_oe 


 %% set the current orbit vector 
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


orbit_current = [sm, ecc, inc, RA, ArgPer];

%% set the master orbit vector 

sm_master = master_oe(1);
ecc_master = master_oe(2);
RA_master = master_oe(3);
inc_master = master_oe(4);
ArgPer_master = master_oe(5);
TrueA_master = master_oe(6);
T_master = master_oe(7);
b_master = master_oe(8);
p_master = master_oe(9);
r_master = master_oe(10);
h_master = master_oe(11);

master_orbit = [sm_master, ecc_master, inc_master, RA_master, ArgPer_master];

%% set the deltas

err = orbit_current - master_orbit 


K = gains (orbit, 1);

%% calculate the K*err 


Kerr(1) = err(1)*K(1);
Kerr(2) = err(2)*K(2);
Kerr(3) = err(3)*K(3);
Kerr(4) = err(4)*K(4);
Kerr(5) = err(5)*K(5);










