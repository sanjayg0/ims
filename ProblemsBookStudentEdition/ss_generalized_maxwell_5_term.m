function ss_generalized_maxwell_5_term()
% Usage: ss_generalized_maxwell_5_term()
% Purpose: Plot steady state response 
%          of 5-term sls Maxwell model

    % Set material parameters
    E0      = 3.35;     % MPa
    E1      = 322.35;   % MPa
    tau1    = 2.3e-3;   % s
    E2      = 129.84;   % MPa
    tau2    = 2.3e-2;   % s
    E3      = 45.12;    % MPa
    tau3    = 0.16;     % s
    E4      = 12.12;    % MPa
    tau4    = 1.12;     % s
    E5      = 5.55;     % MPa
    tau5    = 22.95;    % s

    % Compute storage modulus, loss modulus, 
    % and tan delta as a function of frequency
    omega    = logspace(-1, 1, 100);
    N        = length(omega);
    E_p      = zeros(N,1);
    E_pp     = zeros(N,1);
    tandelta = zeros(N,1);

    for i= 1:N
        E_p(i) = E0 + ...
                E1*(omega(i)*tau1)^2/(1+(omega(i)*tau1)^2) + ...
                E2*(omega(i)*tau2)^2/(1+(omega(i)*tau2)^2) + ...
                E3*(omega(i)*tau3)^2/(1+(omega(i)*tau3)^2) + ...
                E4*(omega(i)*tau4)^2/(1+(omega(i)*tau4)^2) + ...
                E5*(omega(i)*tau5)^2/(1+(omega(i)*tau5)^2) ;

        E_pp(i) = ...
                E1*(omega(i)*tau1)/(1+(omega(i)*tau1)^2) + ...
                E2*(omega(i)*tau2)/(1+(omega(i)*tau2)^2) + ...
                E3*(omega(i)*tau3)/(1+(omega(i)*tau3)^2) + ...
                E4*(omega(i)*tau4)/(1+(omega(i)*tau4)^2) + ...
                E5*(omega(i)*tau5)/(1+(omega(i)*tau5)^2) ; 

        tandelta(i) = E_pp(i)/E_p(i);    
    end    

    % Plot steady state functions
    figure(1);
    yyaxis left;
    semilogx(omega,E_p,'k-',omega,E_pp,'b-','LineWidth',2);
    grid on
    xlabel('$\omega$, rad/s','FontSize',20,...
        'Interpreter','latex'); 
    ylabel('$E^\prime,\  E^{\prime\prime}$\, MPa',...
        'FontSize',20,'Interpreter','latex'); 
    set(gca,'XMinorTick','On', 'Fontsize', 16)
    set(gca,'YMinorTick','On', 'Fontsize', 16)
    set(gca,'Xlim',[.1,10])
    set(gca,'Ylim',[0,100])
    hold on;
    yyaxis right;
    semilogx(omega,tandelta,'r','LineWidth',2);
    ylabel('$\tan \delta$','FontSize',20,'Interpreter','latex'); 
    set(gca,'YMinorTick','On', 'Fontsize', 16)
    set(gca,'Ylim',[0,1.25])
    leg = legend('$E^\prime$','$E^{\prime\prime}$',...
        '$\tan\delta$');
    set(leg,'FontSize',20,'Interpreter','latex',...
        'Location','SouthEast')
    hold off 
end