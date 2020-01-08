function a = AccelThrust(Y,mass,d) 

global AuxParam current_orbit 
%Calculation of DCM for LVLH conversion.

QXx = LVLHTransformation(Y);

%% OLD

% % % thrust force 
% % 
% % thrust = AuxParam.ThrustMag;
% % 
% % %Thrust force magnitude in N
% % %The force is applied in Y axis of LVLH frame, axis similar to velocity.
% % 
% % if d == 0
% %     Thrustlvlh = [ 0 , -thrust , 0];
% % else
% %     Thrustlvlh = [ 0 , thrust , 0];
% % end
% % 
% % %Transformation of Thrust vector to ECI
% % 
% % Thrusteci = inv(QXx) * Thrustlvlh';

%% call Lyapunov for thrust force calculation

T = Lyapunov (current_orbit); 

% Tt = [0, T(2), 0]';

%% Transformation of Thrust vector to ECI

Thrust = inv(QXx) * T;

% Acceleration
a = 1/mass*Thrust;