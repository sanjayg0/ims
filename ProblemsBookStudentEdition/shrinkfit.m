% shrinkfit.m
% Set properties (mm, MPa)
a     = 150;
c     = 200;
b     = 250;
delta = 0.1;
E     = 200.e3;
nu    = 0.28;
p     = 140;
phi   = (1-nu)/(1+nu);

% Evaluate displacement coefficients
fac1 = (b^2-c^2)/(b^2-a^2);
A1   = -fac1*(phi/(1+phi))*delta/c;
B1   = -a^2*fac1*(1/(1+phi))*delta/c;

fac2 = (c^2-a^2)/(b^2-a^2);
A2   = fac2*(phi/(1+phi))*delta/c;
B2   = b^2*fac2*(1/(1+phi))*delta/c;

B3   = p*((1+ nu)/E)*(a^2)*(b^2)/(b^2-a^2);
A3   = (1/b^2)*phi*B3;

% Evaluate stresses in compound cylinder
r1   = linspace(a,c,25);
r1sq = r1.^2;
r2   = linspace(c,b,25);
r2sq = r2.^2;
r3   = linspace(a,b,50);
r3sq = r3.^2;

% Hoop stress in inner/outer/combined cylinder
stt1   = (E/(1+nu))*((A1/phi)+ B1 ./ r1sq );
stt2   = (E/(1+nu))*((A2/phi)+ B2 ./ r2sq );
sttpre = [stt1,stt2];                     

% Hoop stress in monolithic cylinder
stt3   = (E/(1+nu))*((A3/phi)+ B3 ./ r3sq );

% Hoop stress in pressurized compound cylinder
stt    = sttpre + stt3;

% Plot
h = plot(r3,stt,'k',r3,sttpre,'r--',r3,stt3,'b--');
set(h,'LineWidth',2);
axis([150,250,-100,400.]);
axis(['square']);
xlabel('$r$, mm','FontSize',20,'Interpreter','latex')
ylabel('$\sigma_{\theta\theta}$, MPa','FontSize',20,...
	'Interpreter','latex'); 
legend('Final stress distribution',...
	'Shrunk-fit, no internal pressure',...
	'Monolithic, pressurized',...
	'FontSize',16,'Interpreter','latex',...
	'location', 'northeast');
set(gca, 'TickLabelInterpreter', 'latex','XMinorTick','on',...
    'YMinorTick','on', 'Fontsize', 20)