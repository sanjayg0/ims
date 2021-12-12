% Example quasi-isotropic laminate p.364
% T300/5208, h0 = 125um, [0_2, 90_2, 45_2, -45_2]_S

moduli.e1   = 181e9;  % [Pa]
moduli.e2   = 10.3e9;
moduli.nu12 = 0.28;
moduli.g12  = 7.17e9;

layup.h0   = 125e-6;  % [m]
layup.code = [0,0,90,90,45,45,-45,-45,-45,-45,45,45,90,90,0,0];  %[degrees]

laminate(moduli,layup,'qiso_ex.out');