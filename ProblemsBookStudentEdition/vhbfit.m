% vhbfit.m
% Fit experimental stretch-stress data for VHB 4910
  
% Load data
dataset = dlmread('vhb_data.txt');
stretch = dataset(:,1);
s       = dataset(:,2); % kPa

% Stretches for plotting
lambda     = linspace(1,9.5,100);
factor     = @(x) x - 1./(x.*x);

% Fit neo-Hookean model
opts = fitoptions('Method','NonlinearLeastSquares',...
    'Lower',[0],'Upper',[100],'StartPoint',[12]);

neoh = @(a,x) a*factor(x);
f    = fittype(neoh,'options',opts);

[coef,gof] = fit(stretch(1:4),s(1:4),f);
fprintf('neo-Hookean G0 = %e\n\n',coef.a);

% Computer neo-Hookean stresses
stress_neo = neoh(coef.a,lambda);

% Fit Arruda-Boyce model
opts = fitoptions('Method','NonlinearLeastSquares',...
    'Lower',[3,1],'Upper',[20,100],'StartPoint',[10,10]);

lamb = @(x) sqrt((x.*x+2./x)/3);
invl = @(x) x.*(3-x.*x)./(1-x.*x);
ab   = @(a,b,x) b*(a./(3*lamb(x))).*invl(lamb(x)/a).*factor(x);
f    = fittype(ab,'options',opts);

[coef,gof] = fit(stretch,s,f);
fprintf('Arruda-Boyce G0 = %e\n',coef.b);
fprintf('Arruda-Boyce lamL = %e\n\n',coef.a);

% Compute Arruda-Boyce stresses
stress_ab = ab(coef.a,coef.b,lambda);

% Fit Gent model
opts = fitoptions('Method','NonlinearLeastSquares',...
    'Lower',[3,1],'Upper',[300,100],'StartPoint',[250,10]);

gent = @(a,b,x) b*(a./(a-(x.*x+2./x-3))).*factor(x);
f    = fittype(gent,'options',opts);

[coef,gof] = fit(stretch,s,f);

fprintf('Gent G0 = %e\n',coef.b);
fprintf('Gent Im = %e\n',coef.a);

% Compute Gent stresses
stress_gent = gent(coef.a,coef.b,lambda);

% Plot reults
plot(stretch,s,'kd',...
 lambda,stress_neo,'k-.',lambda,stress_ab,'b-',...
 lambda, stress_gent,'r-',...
      'LineWidth',2,'MarkerSize',8,'MarkerFaceColor','k');
grid on
axis([1,10,0,300]);
xlabel('Stretch $\lambda$','FontSize',20,'Interpreter','latex');
ylabel('Engineering stress  s (kPa)',...
    'FontSize',20,'Interpreter','latex');
title('VHB 4910','FontSize',16)
set(gca,'XMinorTick','On', 'Fontsize', 16)
set(gca,'YMinorTick','On', 'Fontsize', 16)
leg = legend('Experiment','Neo-Hookean','Arruda-Boyce','Gent');
set(leg,'FontSize',20,'Interpreter','latex','Location','Best')