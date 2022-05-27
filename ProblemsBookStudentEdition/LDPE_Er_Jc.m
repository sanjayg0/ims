% LDPE_Er_Jc.m
% Script for analyzing stress relaxation and creep in LDPE.

% Import Relaxation data
dataset = dlmread('LDPE_Er.txt');
time    = dataset(1:end,1);      % s
E_r     = dataset(1:end,2)*1E-9; % GPa

% Zero the time  
time  = time - time(1);          % s

% Plot relaxation data
figure(1);
plot(time,E_r,'ko','LineWidth',1);
axis([0, 300,0.1,0.3]); 
xlabel('$t$ (s)','FontSize',16,'Interpreter','Latex');
ylabel('$E_r$ (GPa)','FontSize',16,'Interpreter','Latex');
set(gca, 'TickLabelInterpreter', 'latex',...
    'XMinorTick','on','YMinorTick','off',...
    'Fontsize', 20)
hold on

% Find E_rg 
E_rg = E_r(1) 

% Nonlinear-least-squares fit to
% E_r= E_rg*exp(-(t/tau)^beta)

% Two parameters: a = tau, b=beta
s = fitoptions('Method','NonlinearLeastSquares',...
               'Lower',[1.,0.005],...
               'Upper',[10000.,1.],...
               'Startpoint',[10.,0.5]);
f = fittype('n*exp(-(x/a)^b)','problem','n','options',s);

% Output in c1 with c1.a = tau and c1.b = beta
[c1,gof1] = fit(time,E_r,f,'problem',E_rg);

tau  = c1.a
beta = c1.b

% Compute the Fit
E_r_fit = E_rg*exp(-(time./tau).^beta);

% Add to plot
plot(time,E_r_fit,'r','LineWidth',2)
leg = legend('Experiment','Fit');
set(leg,'FontSize',18,'Location','NorthEast',...
  'Interpreter','Latex')
hold off

% Load creep data
dataset = dlmread('LDPE_Jc.txt');
time    = dataset(1:end,1);      % s
J_c     = dataset(1:end,2)*1.e9; % 1/GPa

% Zero the time
time  = time - time(1);          % s

% Plot relaxation data
figure(2);
plot(time,J_c,'ko','LineWidth',1);

axis([0, 300,0,8]); 
xlabel('t (s)','FontSize',16,'Interpreter','Latex');
ylabel('$J_c$ (GPa$^{-1}$)','FontSize',16,'Interpreter','Latex');
set(gca, 'TickLabelInterpreter', 'latex',...
    'XMinorTick','on','YMinorTick','off',...
    'Fontsize', 20)
hold on

% Estimate J_ce
J_ce = J_c(end); 
 
% Nonlinear-least-squares fit to:
% J_c= J_ce*(1 - exp(-(t/tau)^beta))

% Three params: a = tau, b=beta, c=J_ce
s = fitoptions('Method','NonlinearLeastSquares',...
               'Lower',[0.1,0.005,6],...
               'Upper',[100000.,5.,10],...
               'Startpoint',[10,0.7,J_ce]);
f = fittype(@(a,b,c,x) c*(1-exp(-(x/a).^b)));

% Output in c1 with c1.a = tau, c1.b = alpha, c1.c = J_ce  
[c1,gof1] = fit(time,J_c,f,s);

tau  = c1.a
beta = c1.b
J_ce = c1.c
 
% Compute the Fit
J_c_fit = J_ce.*(1- exp(-(time/tau).^beta));

% Add to plot
plot(time,J_c_fit,'r-','LineWidth',2)
leg = legend('Experiment','Fit');
set(leg,'FontSize',16,'Location','SouthEast',...
  'Interpreter','Latex')
hold off