function T = Lyapunov(orbit)

%% Lyapunov function: calculates the accel-thrust vector 
%% inputs : current orbital elements 



% calculate input matrix G 

G = input_matrix (orbit);

% calculate error 

Kerr = error_function (orbit);

% % calculate gain
% 
% K = gains (orbit);


%% calculate thrust vector 

T = - G' * Kerr';









