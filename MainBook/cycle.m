function [x,t] = cycle(xmax,xmin,xdot,numrev,Nstep)
% Usage: [x,t] = cycle(xmax,xmin,xdot,numrev,Nstep) 
%
% Purpose: returns sawtooth function values and times, starts from zero
%
% Input: xmax   -- max value
%        xmin   -- min value
%        xdot   -- slope of the sawtooth
%        numrev -- number of reversals (no. half-periods)
%        Nstep  -- time steps in a half-cycle (except first quarter-cycle)
%
% Output: x -- vector of values
%         t -- vector of times

    % Error checking
    if nargin ~= 5
       error('Not enough input arguments.');
    elseif xmax == xmin
        error('xmax must not be equal to xmin');
    elseif numrev < 0
        error('numrev must be greater than or equal to 0');
    elseif Nstep < 2
        error('Nstep must be greater than or equal to 2');
    end

    % Fix inputs in case they are not integers
    numrev=floor(numrev);
    Nstep=floor(Nstep);

    % Determine total no. time steps and zero arrays
    N = Nstep * (numrev+1);
    t=zeros(N,1);
    x=zeros(N,1);

    % Initialize quarter-period (t1) and half-period (t2)
    t1 = abs(xmax/xdot);
    t2 = abs((xmax-xmin)/xdot);

    % Set delta-x and delta-t for quarter-cycle and half-cycle
    delt1 = t1/Nstep;
    delx1 = xmax/Nstep;
    delt2 = t2/Nstep;
    delx2 = (xmax-xmin)/Nstep;

    % Compute first quarter-cycle
    n=2;
    for i=1:Nstep
        x(n) = x(n-1) + delx1;
        t(n) = t(n-1) + delt1;
        n=n+1;
    end

    % Compute each reversal
    for i=1:numrev
        % Flip direction
        delx2 = -delx2;

        % Compute values
        for i=1:Nstep
            x(n) = x(n-1) + delx2;
            t(n) = t(n-1) + delt2;
            n=n+1;
        end
    end
end
