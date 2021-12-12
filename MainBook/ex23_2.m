% Example 23.2
% AS4/3501-6, h0 = 125um, [0/90]_S

moduli.e1   = 126e9;  % [Pa]
moduli.e2   = 11e9;
moduli.nu12 = 0.28;
moduli.g12  = 6.6e9;

loads.nx  = -1;  % -1 for compression, +1 for tension, [N/m]
loads.ny  =  0;
loads.nxy =  0;
loads.mx  =  0;  % [N-m/m]
loads.my  =  0;
loads.mxy =  0;

failure.f1t = 1950e6;  % [Pa]
failure.f1c = 1480e6;
failure.f2t = 48e6;
failure.f2c = 200e6;
failure.f12 = 79e6;

layup.h0   = 125e-6;  % [m]
layup.code = [0, 90, 90, 0];  % [degrees]

laminate(moduli,layup,'ex23_2.out',loads,failure);
