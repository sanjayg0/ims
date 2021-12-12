% Example balanced-ply laminate p.363
% T300/5208, h = 125um, [45,-45]_4S

moduli.e1   = 181e9;  % [Pa]
moduli.e2   = 10.3e9;
moduli.nu12 = 0.28;
moduli.g12  = 7.17e9;

layup.h0   = 125e-6;  % [m]
layup.code = [45,-45,45,-45,45,-45,45,-45,-45,45,-45,45,-45,45,-45,45];  % [degrees]

laminate(moduli,layup,'balply_ex.out');