% strainlife1.m
% Set plot range
x=logspace(0,8,100); % 2N_f

% 1015 steel
E         = 206000;  % MPa  
sigfprime = 827;     % MPa
b         = -0.11; 
epsfprime = 0.95;  
c         = -0.64;  

% Plot 1015
y = (sigfprime/E)*x.^b + epsfprime*x.^c; %eps_a
h = loglog(x,y,'k-');
set(h,'LineWidth',2);
hold on

% 1045 steel
E         = 206000;   % MPa  
sigfprime = 2274;     % MPa
b         = -0.08; 
epsfprime = 0.25;  
c         = -0.68;  

% Plot 1045
y = (sigfprime/E)*x.^b + epsfprime*x.^c; %eps_a
h = loglog(x,y,'b-');
set(h,'LineWidth',2);
hold off

axis([1, 1.E8,1.e-4,1]);
xlabel('$2N_f$, reversals to failure',...
    'FontSize',20,'Interpreter','latex')
ylabel('$\epsilon_a$, strain amplitude',...
    'FontSize',20,'Interpreter','latex'); 
legend('$1015$ normalized','$1045$ quenched \& tempered',...
    'FontSize',20,'Interpreter','latex',...
    'location', 'northeast');
set(gca, 'TickLabelInterpreter', 'latex','XMinorTick','on',...
    'YMinorTick','of', 'Fontsize', 20)
