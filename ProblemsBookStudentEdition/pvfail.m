% pvfail.m

% First ply failure pressure vessel 
% AS4/3501-6, h0 = 125um, [+/-(alpha)]_S

moduli.e1   = 126e9;  % [Pa]
moduli.e2   = 11e9;
moduli.nu12 = 0.28;
moduli.g12  = 6.6e9;

loads.nx  =  1/2;  % axial stress [N/m]
loads.ny  =  1;    % hoop stress [N/m]
loads.nxy =  0;
loads.mx  =  0;  % [N-m/m]
loads.my  =  0;
loads.mxy =  0;

failure.f1t = 1950e6;  % [Pa]
failure.f1c = 1480e6;
failure.f2t = 48e6;
failure.f2c = 200e6;
failure.f12 = 79e6;

layup.h0   = 125e-6; % [m]

% Loop over a set of angles
als = linspace(0,90,500);
sf  = zeros(500,1);
i   = 1;
for al = als
  layup.code = [al, -al, -al, al];  % [degrees]
  sf(i)      = laminate_mod(moduli,layup,'pvfail.out',...
      loads,failure);
  i = i + 1;
end

% Plot SF versus winding angle
h = plot(als,sf,'Linewidth',2);
xlabel('Winding angle $\alpha$ (deg)','Interpreter','latex');
ylabel('Tsai-Wu safety factor','Interpreter','latex');
set(gca,'FontSize',20)

% Find best winding angle
[sfm, idx] = max(sf);
fprintf('Max SF at alpha = %4.2f\n',als(idx));
