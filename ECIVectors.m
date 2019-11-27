function [M , N] = ECIVectors


%--------------------------------------------------------------------------
%   INPUT: Vector's file from STK 
%   OUTPUT: Matrix M, which contains the ECI vectors
%--------------------------------------------------------------------------


M = dlmread('ICEYE_X2_Matlab ECI Position Vectors_2 .txt');

N = dlmread('ICEYE_X2_Matlab ECI Velocity Vectors_1 .txt');


% for i = 1:1:1588 
%     Xeci1(i) = M(i,1);
%     Yeci1(i) = M(i,2);
%     Zeci1(i) = M(1,3);
% end
% 
% Xeci = Xeci1';
% Yeci = Yeci1';
% Zeci = Zeci1';
% 
% for j = 1:1:1588 
%     X(j,1) = Xeci(j);
%     X(j,2) = 0;
%     X(j,3) = 0;
%     
%     Y(j,1) = 0;
%     Y(j,2) = Yeci(j);
%     Y(j,3) = 0;
% 
%     Z(j,1) = 0;
%     Z(j,2) = 0;
%     Z(j,3) = Zeci(j);
% end
    
end
