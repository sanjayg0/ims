function [g,dg] = rdresid(X,arg)
% Usage: [g,dg] = rdresid(X,arg)
%
% Purpose: Compute rate dependent plasticity residual function and its
%          derivative wrt strain-rate
%
% Input: X   -- tensile equivalent plastic strain rate
%        arg -- array of material parameters
%
% Output: g  -- value of residual function
%         dg -- derivative of g with respect to equiv. plas. strain-rate

    % Extract material parameters
    sigtr = arg(1);
    Yn    = arg(2);
    Hn    = arg(3);
    delt  = arg(4);
    E     = arg(5);
    edot0 = arg(6);
    mrate = arg(7);
    Y0    = arg(8);
    H0    = arg(9);
    Ys    = arg(10);
    rhard = arg(11);

    % Evaluate residual and derivative
    g  = abs(sigtr) -E*delt*X - (Yn + Hn*delt*X)*(X/edot0)^mrate;
    dg = -E*delt - Hn*delt*(X/edot0)^mrate ...
        - (mrate/edot0)*(Yn + Hn*delt*X)*(X/edot0)^(mrate-1);
end
