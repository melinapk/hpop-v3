function reduce_flag = calculate_thrust(mass, delta_sm , mean_sm , delta_t )

global AuxParam

reduce_flag = 0;
mu = 398600.4418e9; 

% calculate the thrust force ( Fe ) with semi-major axis ( sm )

F_thrust = mass*sqrt(mu/power(mean_sm,3))/(2*delta_t)*delta_sm;

%Indicate the reduction of Thrust force

Delta_Thrust = AuxParam.ThrustMag - F_thrust; 

% Change the global parameter

AuxParam.ThrustMag = F_thrust; 

if Delta_Thrust  < 0 
   reduce_flag = 1;
end