% Example 23.1
% T300/5208, h = 125um, [0, 30, -30, 90]

moduli.e1   = 181e9;  % [Pa]
moduli.e2   = 10.3e9;
moduli.nu12 = 0.28;
moduli.g12  = 7.17e9;

layup.h0   = 125e-6;  % [m]
layup.code = [0, +30, -30, 90];  % [degrees]

laminate(moduli,layup,'ex23_1.out');