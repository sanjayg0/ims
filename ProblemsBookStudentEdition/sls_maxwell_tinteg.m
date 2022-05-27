function sls_maxwell_tinteg()
% Usage: sls_maxwell_tinteg()
% Purpose: Driver function for time integration 
%          of Standard-Linear-Solid Maxwell model
         
    % Select loading case: 1 saw, 2 sine, or 3 ramp
    pt = 3;
    
    % Set material parameters
    E0      = 1; % MPa
    E1      = 1; % MPa
    tau1    = 1; % s
    const1  = [E0 E1 tau1];

    switch pt
        case 1 
            % saw-tooth strain input
            edot          = 0.01; % strain rate
            emax          = 0.01; % max strain
            num_reversals = 6;    % num of  reversals
            Nstep         = 10;   % num of incrs in a half-cycle

            % Generate saw-tooth strain history
            [e,t]   = cycle(emax,-emax,edot,num_reversals,Nstep);
            % Create output
            [sigma] = main_sls(const1,e,t); 

        case 2 
            % sinusoidal strain input
            ea       = 0.01;               % strain amplitude
            T_per    = 4;                  % time period
            omega    = 2*pi/T_per;         % angular frequency
            t        = linspace(0,10,100); % time

            % Generate sinusoidal strain history
            e       = ea*sin(omega*t); 
            % Create output
            [sigma] = main_sls(const1,e,t);             

        case 3 
            % ramp strain input
            ea       = 0.01;                      % strain amp 
            t0       = 0.0;                       % ramp time   
            t1       = 1.0;                       % ramp time
            t2       = 9.0;                       % ramp time  
            t3       = 10.0;                      % ramp time
            time1    = linspace(t0,t1,11);        % time1
            time2    = linspace(t1,t2,40);        % time2
            time3    = linspace(t2,t3,10);        % time2
            t        = [time1,time2,time3];       % time

            % Generate ramp strain history
            strain1  = ea*(time1-t0)/(t1-t0);     % strain1
            Ntime    = length(time2);             % len time2
            strain2  = ea*ones(1,Ntime);          % strain2
            strain3  = ea*(t3-time3)/(t3-t2);     % strain3
            e        = [strain1,strain2,strain3]; % strain
            % Create output
            [sigma] = main_sls(const1,e,t);                                  
    end

    % Plot strain-time curve
    figure(1);
    plot(t,e,'k-o','LineWidth',2);
    xlabel('$t$ (s)','FontSize',18,'Interpreter','latex');
    ylabel('$\epsilon$','FontSize',18,'Interpreter','latex');
    set(gca,'XMinorTick','On', 'Fontsize', 16)
    set(gca,'YMinorTick','On', 'Fontsize', 16)

    % Plot stress-time curve
    figure(2);
    plot(t,sigma,'r-o','LineWidth',2);
    xlabel('$t$ (s)','FontSize',18,'Interpreter','latex');
    ylabel('$\sigma$ (MPa)','FontSize',18,'Interpreter','latex');
    set(gca,'XMinorTick','On', 'Fontsize', 16)
    set(gca,'YMinorTick','On', 'Fontsize', 16)

    % Plot strain-time and stress-time curves together
    figure(3);
    yyaxis left;
    plot(t,e,'k-o','LineWidth',2);
    xlabel('$t$ (s)','FontSize',20,'Interpreter','latex'); 
    ylabel('$\epsilon$','FontSize',20,'Interpreter','latex'); 
    set(gca,'XMinorTick','On', 'Fontsize', 16)
    set(gca,'YMinorTick','On', 'Fontsize', 16)
    hold on;
    yyaxis right;
    plot(t,sigma,'r-o','LineWidth',2);
    ylabel('$\sigma$ (MPa)','FontSize',20,'Interpreter','latex'); 
    set(gca,'XMinorTick','On', 'Fontsize', 16)
    set(gca,'YMinorTick','On', 'Fontsize', 16)
    legend_handle= legend('$\epsilon$','$\sigma$');
    set(legend_handle,'FontSize',20,'Interpreter','latex',...
        'Location','Best')
    hold off 

    % Plot the stress-strain curve
    figure(4);
    plot(e,sigma,'b-o','LineWidth',2);
    xlabel('$\epsilon$','FontSize',18,'Interpreter','latex');
    ylabel('$\sigma$ (MPa)','FontSize',18,'Interpreter','latex');
    set(gca,'XMinorTick','On', 'Fontsize', 16)
    set(gca,'YMinorTick','On', 'Fontsize', 16)
end

function [sigma] = main_sls(const1,strain,time)
% Usage: [sigma] = main_sls(const1,strain,time)
% Purpose: Time integration procedure given strain history
%
% Input: const1 -- Material property array
%        strain -- Strain history
%        time   -- Time history
%
% Output: sigma -- stress history

    % Extract the material parameters
    E0    = const1(1);
    E1    = const1(2);
    tau1  = const1(3);
    eta1  = tau1*E1;
    
    % Create vectors to store results
    N     = length(strain); % number of total steps
    ev    = zeros(N,1);     % viscous strain   
    sigma = zeros(N,1);     % stress 

    % Perform the incremental time integration
    for n=1:(N-1)
        delt = time(n+1) - time(n);
        dele = strain(n+1) - strain(n);
        %
        % Compute trial elastic strain
        epsetr =  strain(n+1) - ev(n);
        %
        % Update 
        sigma(n+1) = E0*strain(n+1) + (eta1/(tau1+delt))*epsetr; 
        ev(n+1)    = ev(n) + (delt/(tau1+delt))*epsetr;      
    end
end
 
function [x,t] = cycle(xmax,xmin,xdot,numrev,Nstep)
% Usage: [x,t] = cycle(xmax,xmin,xdot,numrev,Nstep)
% Purpose: Returns sawtooth function values and times, 
%          starts from zero
%
% Input: xmax   -- max value
%        xmin   -- min value
%        xdot   -- slope of the sawtooth
%        numrev -- number of reversals (no. half-periods)
%        Nstep  -- time steps in a half-cycle (except first
%                  quarter-cycle)
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
    numrev = floor(numrev);
    Nstep  = floor(Nstep);

    % Determine total no. time steps and zero arrays
    N = Nstep*(numrev+1);
    t = zeros(N,1);
    x = zeros(N,1);

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