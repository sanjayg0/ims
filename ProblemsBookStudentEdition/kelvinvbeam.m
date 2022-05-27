% kelvinvbeam.m
% Set parameters
L   = 1;      % length of beam, m
b   = 0.04;   % edge length of cross-section of the beam, m
E   = 0.5E9;  % Young's modulus Pa
eta = 5.E9;   % viscosity Pa.s
p0  = 0.35E3; % distributed loading N/m.

% Compute factors and set up time-stepping
Imom     = b^4/12;
fac      = (5*p0*L^4)/(384*E*Imom);
tot_time = 200;        % total time
N        = 5000;       % number of increments
dt       = tot_time/N; % increment in time
time     = zeros(N,1); % time
delta    = zeros(N,1); % strain 

% Loop over time to compute displacement
for i=2:N
    time(i) = time(i-1)+dt;
    if time(i) <= 100 
       delta(i)  =  -fac*(1-exp(-time(i)/(eta/E)))*100;
    elseif (time(i)> 100) 
       delta(i)  =  -0.0427*(exp(eta/E)-1)...
                  * exp(-time(i)/(eta/E))*100;
    end 
end 

% Plot result 
plot(time,delta,'k-','LineWidth',2);
axis([0, 200, -5, 0.5]);                  
xlabel('$t$ (s)','FontSize',20,'Interpreter','latex');
ylabel('$\delta$ (cm)','FontSize',20,'Interpreter','latex'); 
set(gca, 'TickLabelInterpreter', 'latex','XMinorTick','on',...
    'YMinorTick','off', 'Fontsize', 20);
