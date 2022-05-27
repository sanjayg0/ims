% rdtinteg1.m
% Script for rate dependent plasticity
% explicit time integration linear loading

% Select case: 1 or 2 
pt = 2;

% Set material parameters
switch pt
    case 1
        E      = 200e9;
        edot0  = 1.E-3;
        mrate  = 0.02;
        Y0     = 250e6;
        Khard  = 1182e6;
        nhard  = 0.47;     
    case 2
        E     = 10e9;
        edot0 = 1.E-3;
        mrate = 0.2;
        Y0    = 12e6;
        Khard = 5e6;
        nhard = 0.8; 
end

% Set strain history parameters
edots         = [0.01, 0.1, 1.0]; % strain rates
emax          = 0.25;             % max strain
Nstep         = 10000;            % Use a large number since 
                                  % FE is an explicit procedure
                                  
% Initialize figure
h      = figure(pt);
ColOrd = get(gca,'ColorOrder');
hold on;

cnt = 1;
for edot = edots    
    % Generate strain history and times
    t = linspace(0,emax/edot,Nstep);
    e = edot*t;

    % Create vectors to store results
    N     = length(e);   % total number of steps
    ep    = zeros(N,1);  % plastic strain  
    ebarp = zeros(N,1);  % equiv. tensile plastic strain 
    Y     = zeros(N,1);  % flow resistance  
    sigma = zeros(N,1);  % stress 

    Y(1)  = Y0;          % initialize Y

    % Perform the explicit time integration
    for n=1:(N-1)      
        delt = t(n+1) - t(n);  

        % Compute the direction of plastic flow
        np(n) = sign(sigma(n));

        % Compute the flow resistance
        Y(n) = Y0 + Khard*(ebarp(n))^nhard;

        % Compute the equiv. tensile plastic strain
        ebarpdot = edot0*(abs(sigma(n))/Y(n))^(1/mrate);        
        debarp   = delt*ebarpdot;

        % Update model values
        ep(n+1)    = ep(n) + debarp*np(n);
        ebarp(n+1) = ebarp(n) + ebarpdot*delt;      
        Y(n+1)     = Y0 + Khard*(ebarp(n+1))^nhard;
        sigma(n+1) = E*(e(n+1)-ep(n+1));    
    end 
    
    % Plot stress-strain curve
    plot(e,sigma/1.E6,'Color',ColOrd(cnt,:),'LineWidth',2);
    cnt = cnt + 1;
end

% Ajust plot
switch pt
    case 1
        axis([0,.25,0,1200]);
    case 2
        axis([0,0.25,0,60]);
end
xlabel('$\epsilon$','Interpreter','latex');
ylabel('$\sigma$ (MPa)','Interpreter','latex');
set(gca,'XMinorTick','on','YMinorTick','on')
set(gca, 'TickLabelInterpreter','latex',...
    'XMinorTick','on','YMinorTick','on','Fontsize',18)

hl = legend('$\dot\epsilon=0.01$/s',...
            '$\dot\epsilon=0.1$/s',...
            '$\dot\epsilon=1.0$/s');
set(hl,'FontSize',16,'Interpreter','latex','Location','Best')

