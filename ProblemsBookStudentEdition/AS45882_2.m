% AS45882_2.m

% Set laminae properties AS4/5882 carbon/epoxy
moduli.e1   = 127.4e9;
moduli.e2   = 10.2e9;
moduli.nu12 = 0.275;
moduli.g12  = 5.59e9;

% Set layup
layup.h0   = 125.0e-6;

% Compute composite properties (i)
layup.code = [0, 0, 45, -45, 45, -45, 90, 90, 90 ,90, -45,...
    45, -45, 45, 0, 0];
laminate(moduli,layup,'p78i.out');

% Compute composite properties (ii)
layup.code = [0, 45, -45, 90, 0, 45, -45, 90, 90, -45, 45,...
    0, 90, -45, 45, 0];
laminate(moduli,layup,'p78ii.out');

% Compute composite properties (iii)
layup.code = [45, -45, 0, 0, 45, -45, 90, 90, 90, 90, -45,...
    45, 0, 0, -45, 45];
laminate(moduli,layup,'p78iii.out');

% Compute composite properties (iv)
layup.code = [45, -45, 0, 90, 45, -45, 0, 90, 90, 0, -45,...
    45, 90, 0, -45, 45];
laminate(moduli,layup,'p78iv.out');

% Compute composite properties (v)
layup.code = [45, -45, 45, -45, 0, 0, 90, 90, 90, 90, 0,...
    0, -45, 45, -45, 45];
laminate(moduli,layup,'p78v.out');

% Compute composite properties (vi)
layup.code = [45, -45, 45, -45, 90, 90, 0, 0, 0, 0, 90,...
    90, -45, 45, -45, 45];
laminate(moduli,layup,'p78vi.out');

% Compute composite properties (vii)
layup.code = [90, 90, 45, -45, 45, -45, 0, 0, 0, 0, -45,...
    45, -45, 45, 90, 90];
laminate(moduli,layup,'p78vii.out');
