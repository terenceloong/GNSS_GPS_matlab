function [K1, K2] = orderTwoLoopCoef(LBW, zeta, k)
% 二阶环路系数
%   Inputs:
%       LBW           - Loop noise bandwidth
%       zeta          - Damping ratio
%       k             - Loop gain
%
%   Outputs:
%       K1, K2        - Loop filter coefficients 

% Solve natural frequency
Wn = LBW*8*zeta / (4*zeta^2 + 1);

% solve for K1 & K2
K1 = 2*zeta*Wn / k;
K2 = Wn^2 / k;

end