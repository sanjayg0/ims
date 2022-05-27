% visco_inv.m
% Script for integrating Jc(t-tau)dEr(tau)/dtau from 0- to t
% from experimental data for Jc and Er.

clear all
close all

% Import data
Er  = importdata('LDPE_Er.txt');
Jc  = importdata('LDPE_Jc.txt');

dt  = mean(diff(Er(:,1))); % Time step between data points

% Calculate numerical derivative of Er at midpoints
dErdt   = diff(Er(:,2))/dt; 

% Create midpoint values of Jc data
Jc_mid  = (Jc(1:end-1,:)+Jc(2:end,:))/2; 

% Define midpoint times
t_mid   = Jc_mid(:,1);  

% Define max value of upper intergation limit
width   = min(length(Jc_mid),length(dErdt));

% Calculate integral numerically
for t=1:1:width
    for tau=1:t
    Jc_back(tau)    = Jc_mid(t-tau+1,2); 
    end
    int(t)  = dot(Jc_back(1:t),dErdt(1:t)*dt);
end


plot(t_mid(1:width), int(1:width),'ko','LineWidth',1)
hold on
plot([-10:1:325],[-10:1:325]>0,'r','LineWidth',2)
axis([-10, 325,0,1.5]); 
xlabel('$t$ (s)','FontSize',16,'Interpreter','Latex');
set(gca, 'TickLabelInterpreter', 'latex',...
    'XMinorTick','on','YMinorTick','off',...
    'Fontsize', 20)
leg = legend('Integral','$h(t)$','Interpreter','Latex');
set(leg,'FontSize',16,'Location','SouthEast',...
  'Interpreter','Latex')
