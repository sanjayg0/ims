% rosette1.m
% Strain gauge rosette program

% Set orientations of the strain-gauge rosettes in degress
theta1 = 0;   
theta2 = 45;
theta3 = 90;

% Set orientations of the strain-gauge rosettes in radians
theta1 = theta1*pi/180;   
theta2 = theta2*pi/180;
theta3 = theta3*pi/180;

% Set experimentally-measured strains
% in the three orientations in microstrain
eps11star = 100;
eps22star = 200;
eps33star = 900;

epsstar   = [eps11star eps22star eps33star]';

% Construct the A-matrix
A = [cos(theta1)^2 sin(2*theta1) sin(theta1)^2;
     cos(theta2)^2 sin(2*theta2) sin(theta2)^2;
     cos(theta3)^2 sin(2*theta3) sin(theta3)^2];

% Calculate the strain components
% (\epsilon_{11},\epsilon_{12},\epsilon_{22})
eps = A\epsstar;

% Display Results
fprintf('epsilon11 = %5.2f\n',eps(1)); 
fprintf('epsilon12 = %5.2f\n',eps(2)); 
fprintf('epsilon22 = %5.2f\n',eps(3)); 
