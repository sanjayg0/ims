% strainlife2.m
% Select material
material = 1045;
fprintf('Material %d\n',material);

switch (material)
    case 1045
        E         = 206000;  % MPa
        sigfprime = 2274;    % MPa
        b         = -0.08; 
        epsfprime = 0.25;  
        c         = -0.68; 
    case 1015
        E         = 206000;  % MPa
        sigfprime = 827;     % MPa
        b         = -0.11; 
        epsfprime = 0.95;  
        c         = -0.64;
end

% Set mean stress
sig_m = 0; % MPa

for eps_a = [0.002, 0.005, 0.01]
    % Let x = 2N_f --- the number of reversal
    % Define the implicit function to be solved:
    F = @(x) ((sigfprime-sig_m)/E)*(x^b)+epsfprime*(x^c)-eps_a;
    
    % Use the MATLAB function fzero to obtain the solution
    % Restrict the initial guess to the interval [1 1.E10]
    reversals = fzero(F,[1.0 1.E10]); % number of reversals
    fprintf('Strain amp = %e\t reversals = %e\n',...
        eps_a, reversals);
end
