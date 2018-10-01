function inside = D(x) 
%--------------------------------------------------------------------------
% Matlab M-file Project: HyEQ Toolbox @  Hybrid Systems Laboratory (HSL), 
% https://hybrid.soe.ucsc.edu/software
% http://hybridsimulator.wordpress.com/
% Filename: D.m
%--------------------------------------------------------------------------
% Description: Jump set
% Return 0 if outside of D, and 1 if inside D
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%   See also HYEQSOLVER, PLOTARC, PLOTARC3, PLOTFLOWS, PLOTHARC,
%   PLOTHARCCOLOR, PLOTHARCCOLOR3D, PLOTHYBRIDARC, PLOTJUMPS.
%   Copyright @ Hybrid Systems Laboratory (HSL),
%   Revision: 0.0.0.3 Date: 05/20/2015 3:42:00

global cTilde_10 d_10 gamma_0 r0Prime

% state
z1 = x(1);
z2 = x(2);
q = x(3);

normGradL = abs(GradientL(z1));
halfZ2Squared = (1/2)*z2^2;

V0 = gamma_0 * alpha1(distance(z1)) + (1/2)*z2^2;

% Determine whether system is inside the jump set or not
if (q == 0 && (V0 >= r0Prime))||(q == 1 && (normGradL <= cTilde_10 && halfZ2Squared <= d_10))  
    inside = 1;
else
    inside = 0;
end

end