function ritinteg()
% Usage: ritinteg
%
% Purpose: Driver function for rate independent plasticity
%          time integration

    % Set material parameters
    E      = 200e9;
    Y0     = 250e6;
    H0     = 2000e6;
    Ys     = 500e6;
    r      = 1.0;
    const1 = [E Y0 H0 Ys r];

    % Set strain history parameters
    edot          = 0.001; % strain rate
    emax          = 0.02;  % max strain
    num_reversals = 20;    % number of load reversals
    const2        = [edot emax num_reversals];

    % Call function to compute stress-strain history
    [e, sigma, t] = main_ri(const1, const2);

    % Plot the stress-strain curve
    figure;
    h = plot(e,sigma/1.E6,'k');
    set(h,'LineWidth',1);
    axis([-0.025,0.025,-600,600]);
    xlabel('$\epsilon$','Interpreter','latex');
    ylabel('$\sigma$, MPa','Interpreter','latex');
    set(gca,'FontSize',16)
    set(gca,'XMinorTick','On')
    set(gca,'YMinorTick','On')
end
