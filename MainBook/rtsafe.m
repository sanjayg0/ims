function [x] = rtsafe(x1,x2,xacc,arg,maxit,fhandle)
% Usage: [x] = rtsafe(x1,x2,xacc,arg,maxit,fhandle)
%
% Purpose: Compute the root of a non-linear function using an interval
%          Newton method
%
% Input: x1      -- lower bound to interval with root
%        x2      -- upper bound to interval with root
%        xacc    -- interval tolerance for root
%        arg     -- array of arguments for non-linear function
%        maxit   -- maximum number of iterations
%        fhandle -- function handle to non-linear equation
%
% Output: x -- root to non-linear equation

    % Compute function values at upper and lower bounds
    [fl,df] = fhandle(x1,arg);
    [fh,df] = fhandle(x2,arg);

    % Check end point of interval for root
    if (fl == 0)
        x = x1;
        return;
    elseif (fh == 0)
        x = x2;
        return;
    end

    % Determine end point of interval with negative value
    if (fl*fh > 0)
        error('rtsafe: interval does not necessarily bound root');
    elseif (fl < 0)
        xl = x1;
        xh = x2;
    else
        xh = x1;
        xl = x2;
    end

    % Compute midpoint and interval size
    x     = 0.5 * (x1 + x2);
    dxold = abs(x2 - x1);
    dx    = dxold;

    % Evaluate at midpoint
    [f,df] = fhandle(x,arg);

    % Iterate to find root
    for j=1:maxit
        % Compute tests for bounding properties of Newton step 
        test1 = ((x - xh)*df - f) * ((x - xl)*df - f);
        test2 = abs(2*f) - abs(dxold * df);

        if (test1 > 0 || test2 > 0)
            % Compress interval for unsafe Newton step
            dxold = dx;
            dx    = 0.5 * (xh - xl);
            x     = xl + dx;
            if (xl == x)
                return;
            end
        else
            % Newton step if safe
            dxold = dx;
            dx    = f/df;
            temp  = x;
            x     = x - dx;
            if (temp == x)
                return;
            end
        end

        % Check interval tolerance
        if (abs(dx) < xacc)
            return;
        end

        % Get value for next step
        [f,df] = fhandle(x,arg);

        % Determine side containing the root
        if (f < 0)
            xl = x;
        else
            xh = x;
        end
    end
    error('rtsafe: exceeded maximum iterations');
end
