function generalized_maxwell_5_term_tinteg()
% Usage: generalized_maxwell_5_term_tinteg()
% Purpose: Driver function for time integration 
%          of 5-term sls Maxwell model

    % Select loading case 1, 2, or 3
    pt = 1;

    % Set material parameters
    E0      = 3.35;     % MPa
    E1      = 322.35;   % MPa
    tau1    = 2.3E-3;   % s
    E2      = 129.84;   % MPa
    tau2    = 2.3E-2;   % s
    E3      = 45.12;    % MPa
    tau3    = 0.16;     % s
    E4      = 12.12;    % MPa
    tau4    = 1.12;     % s
    E5      = 5.55;     % MPa
    tau5    = 22.95;    % s
    const   = [E0 E1 tau1 E2 tau2 E3 tau3 E4 tau4 E5 tau5];

    switch pt
        % Stress-strain curves at different strain rates
        case 1  
             % first strain rate
             edot          = 0.01; % strain rate
             emax          = 0.02; % max strain
             num_reversals = 1;    % num of  reversals
             Nstep         = 100;  % num of incs in a half-cycle

             % Generate strain history
             [e,t] = cycle(emax,-emax,edot,num_reversals,Nstep);

             % Create output
             [sigma] = main_gen_max(const,e,t);             

             % Plot the stress-strain curve
             figure(1)
             plot(e,sigma,'k-','LineWidth',2);
             axis([0,0.025, -1.,1.5]);        
             xlabel('$\epsilon$','FontSize',18,...
                 'Interpreter','latex');
             ylabel('$\sigma$, MPa','FontSize',18,...
                 'Interpreter','latex');
             set(gca,'XMinorTick','On', 'Fontsize', 16);
             set(gca,'YMinorTick','On', 'Fontsize', 16);
             grid on
             hold on

             % second strain rate
             edot          = 0.05; % strain rate

             % Generate saw-tooth strain history
             [e,t] = cycle(emax,-emax,edot,num_reversals,Nstep);

             % Create output
             [sigma] = main_gen_max(const,e,t);             

             % Plot the stress-strain curve  
             plot(e,sigma,'r-','LineWidth',2);

             % third strain rate            
             edot          = 0.1; % strain rate

             % Generate saw-tooth strain history
             [e,t] = cycle(emax,-emax,edot,num_reversals,Nstep);

             % Create output
             [sigma] = main_gen_max(const,e,t);             

             % Plot the stress-strain curve  
             plot(e,sigma,'b-','LineWidth',2);
             legend_handle= legend(...
                '$\dot\epsilon=0.01\, \mathrm{s}^{-1}$ ',...
                '$\dot\epsilon=0.05\, \mathrm{s}^{-1}$ ',...
                '$\dot\epsilon=0.1\, \mathrm{s}^{-1}$ ');
             set(legend_handle,'FontSize',20,'Interpreter',...
                'latex','Location','Best');
             hold off
             
        % Stress-strain curves at different strain levels
        case 2  
             % first strain level
             edot          = 0.1;  % strain rate
             emax          = 0.01; % max strain
             num_reversals = 1;    % num of reversals
             Nstep         = 100;  % num of incs in a half-cycle

             % Generate saw-tooth strain history
             [e,t] = cycle(emax,-emax,edot,num_reversals,Nstep);

             % Create output
             [sigma] = main_gen_max(const,e,t);             

             % Plot the stress-strain curve
             figure(2)
             plot(e,sigma,'k-','LineWidth',2);
             axis([0,0.025, -1.,1.5]);        
             xlabel('$\epsilon$','FontSize',18,...
                 'Interpreter','latex');
             ylabel('$\sigma$, MPa','FontSize',18,...
                 'Interpreter','latex');
             set(gca,'XMinorTick','On', 'Fontsize', 16)
             set(gca,'YMinorTick','On', 'Fontsize', 16)
             grid on
             hold on

             % second strain level
             emax          = 0.015;% max strain

             % Generate saw-tooth strain history
             [e,t] = cycle(emax,-emax,edot,num_reversals,Nstep);

             % Create output
             [sigma] = main_gen_max(const,e,t);             

             % Plot the stress-strain curve  
             plot(e,sigma,'r-','LineWidth',2);

             % third strain level        
             emax          = 0.02;% max strain

             % Generate saw-tooth strain history
             [e,t] = cycle(emax,-emax,edot,num_reversals,Nstep);

             % Create output
             [sigma] = main_gen_max(const,e,t);             

             % Plot the stress-strain curve  
             plot(e,sigma,'b-','LineWidth',2);
             legend_handle= legend(...
                '$\epsilon=0.01$ ',...
                '$\epsilon=0.015$ ',...
                '$\epsilon=0.02$ ');
             set(legend_handle,'FontSize',20,'Interpreter',...
                'latex','Location','Best');
             hold off
             
        % sinusoidal strain input
        case 3 
             % first frequency             
             ea    = 0.02;                   % strain amplitude
             omega = 5;                      % ang. freq. rad/s
             T_per = 2*pi/omega;             % time period
             t     = linspace(0,2*T_per,200);% time
             e     = ea*sin(omega*t);        % strain 

             % Create output
             [sigma] = main_gen_max(const,e,t);   

             % Plot the stress-strain curve
             figure(3);
             plot(e,sigma,'k-','LineWidth',2);
             axis([-0.025,0.025, -3,3]);  
             xlabel('$\epsilon$','FontSize',20,...
                 'Interpreter','latex');
             ylabel('$\sigma$, MPa','FontSize',20,...
                 'Interpreter','latex');
             set(gca,'XMinorTick','On', 'Fontsize', 16)
             set(gca,'YMinorTick','On', 'Fontsize', 16)
             grid on            
             hold on
             
             % second frequency            
             omega = 15;                     % ang. freq. rad/s
             T_per = 2*pi/omega;             % time period
             t     = linspace(0,2*T_per,200);% time
             e     = ea*sin(omega*t);        % strain 

             % Create output
             [sigma] = main_gen_max(const,e,t);   

             % Plot the stress-strain curve  
             plot(e,sigma,'r-','LineWidth',2);

             % third frequency           
             omega = 30;                     % ang. freq. rad/s
             T_per = 2*pi/omega;             % time period
             t     = linspace(0,2*T_per,200);% time
             e     = ea*sin(omega*t);        % strain 
             
             % Create output
             [sigma] = main_gen_max(const,e,t);            
             
             % Plot the stress-strain curve  
             plot(e,sigma,'b-','LineWidth',2);

             legend_handle= legend(...
                '$\omega=1\, \mathrm{rad/s}$ ',...
                '$\omega=10\, \mathrm{rad/s}$ ',...
                '$\omega=100\, \mathrm{rad/s}$ ');
             set(legend_handle,'FontSize',20,'Interpreter',...
                'latex','Location','Best');
             hold off
    end
end

function [sigma] = main_gen_max(const,e,t)
% Usage: [sigma] = main_gen_max(const,strain,time)
% Purpose: finite difference time integration procedure 
%          for an input strain history
%
% Input: const -- Material property array
%
% Output: sigma -- stress history


    % Extract the material parameters
    E0    = const(1);
    E1    = const(2);
    tau1  = const(3);
    E2    = const(4);
    tau2  = const(5);
    E3    = const(6);
    tau3  = const(7);   
    E4    = const(8);
    tau4  = const(9);  
    E5    = const(10);
    tau5  = const(11);   
    eta1  = E1*tau1;
    eta2  = E2*tau2;
    eta3  = E3*tau3;
    eta4  = E4*tau4;
    eta5  = E5*tau5;

    % Create vectors to store results
    N     = length(e);   % total number of steps for all cycles
    ev1   = zeros(N,1);  % viscous strain
    ev2   = zeros(N,1);  % viscous strain 
    ev3   = zeros(N,1);  % viscous strain 
    ev4   = zeros(N,1);  % viscous strain 
    ev5   = zeros(N,1);  % viscous strain     
    sigma = zeros(N,1);  % stress 

    % Perform the incremental time integration
    for n=1:(N-1)
        delt = t(n+1) - t(n);
        dele = e(n+1) - e(n);
        
        % Compute trial elastic strains
        eetr1 =  e(n+1) - ev1(n);
        eetr2 =  e(n+1) - ev2(n);
        eetr3 =  e(n+1) - ev3(n);
        eetr4 =  e(n+1) - ev4(n);
        eetr5 =  e(n+1) - ev5(n);

        % Update stress and viscous strains
        sigma(n+1) = E0*e(n+1) ...
            + (eta1/(tau1+delt))*eetr1 ...
            + (eta2/(tau2+delt))*eetr2 ...
            + (eta3/(tau3+delt))*eetr3 ...
            + (eta4/(tau4+delt))*eetr4 ...
            + (eta5/(tau5+delt))*eetr5;
        ev1(n+1)   = ev1(n) + (delt/(tau1+delt))*eetr1; 
        ev2(n+1)   = ev2(n) + (delt/(tau2+delt))*eetr2; 
        ev3(n+1)   = ev3(n) + (delt/(tau3+delt))*eetr3; 
        ev4(n+1)   = ev4(n) + (delt/(tau4+delt))*eetr4; 
        ev5(n+1)   = ev5(n) + (delt/(tau5+delt))*eetr5;
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
 
