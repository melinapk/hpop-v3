function orbit_flag = Altitude_Control(sm,cta,pta,k)

global AuxParam 

%Control method 1 
%Timing control in order to reach altitude increase of 50Km 

% direction = 0;
% if t > 2236
%    AuxParam.Thrust = 0;
% end
% 
% if t > 84164 
%     AuxParam.Thrust = 1;
%     direction = 1;
% end


%Control method 2 
%Semi - major axis control in order to reach increase of 50Km. 
%Elements used: Semi-major axis (sm) , True Anomaly (TrueA)

orbit_flag = 0;

if sm > 7095
    AuxParam.Thrust = 0;
end
  
if cta < pta && AuxParam.Thrust == 0
   orbit_flag = 1;
end
 
 
 if k > 13
     AuxParam.Thrust = 1; 
     AuxParam.direction = 0; 
 end
 end 
