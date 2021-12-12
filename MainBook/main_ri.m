function [e, sigma, t] = main_ri(const1, const2)
% Usage: [e, sigma, t] = main_ri(const1, const2)
%
% Purpose: Computute rate independent plasticity time histories
%
% Input: const1 -- Material property array
%        const2 -- Strain history parameters
%
% Output: e     -- strain time history
%         sigma -- stress history
%         t     -- array of times

    % Extract material parameters
    E     = const1(1);
    Y0    = const1(2);
    H0    = const1(3);
    Ys    = const1(4);
    rhard = const1(5);

    % Extract strain history parameters
    edot          = const2(1); % strain rate
    emax          = const2(2); % max strain
    num_reversals = const2(3); % number of reversals

    % Specify the number of increments of strain in half-cycle
    Nstep = 1000;

    % Generate sawtooth strain history
    [e,t] = cycle(emax,-emax,edot,num_reversals,Nstep);

    % Create vectors to store results
    N     = length(e);  % total number of steps for all cycles
    ep    = zeros(N,1); % plastic strain  
    Y     = zeros(N,1); % flow resistance  
    sigma = zeros(N,1); % stress 

    Y(1)  = Y0;         % initialize Y

    % Perform the incremental time integration
    for n=1:(N-1)
        delt  = t(n+1) - t(n);
        dele  = e(n+1) - e(n);
        sigtr = E * (e(n+1) - ep(n));
        ftr   = abs(sigtr) - Y(n);

        if ftr <= 0
            % Elastic step
            sigma(n+1) = sigtr;        
            ep(n+1)    = ep(n);
            Y(n+1)     = Y(n);
        else
            % Plastic step integration with semi-implicit method
            Hn       = H0*(1-Y(n)/Ys)^rhard;
            delepbar = (abs(sigtr)-Y(n))/(E+Hn);

            % Update
            ep(n+1)    = ep(n) + delepbar*sign(sigtr);
            Y(n+1)     = Y(n)  + Hn*delepbar;
            sigma(n+1) = sigtr - E * delepbar*sign(sigtr);
        end
    end
end
