function rdtinteg()
% Usage: rdtinteg
%
% Purpose: Driver function for rate dependent plasticity
%          time integration

    % Select case: 1 or 2 for Example 10.3
    pt = 2;

    % Set material parameters
    switch pt
        case 1
            E      = 200e9;
            edot0  = 1.E-3;
            mrate  = 0.02;
            Y0     = 250e6;
            H0     = 2000e6;
            Ys     = 500e6;
            rhard  = 1.0;     
        case 2
            E     = 20e9;
            edot0 = 1.E-3;
            mrate = 0.2;
            Y0    = 10e6;
            H0    = 100e6;
            Ys    = 20e6;
            rhard = 1.0;
    end
    const1 = [E edot0 mrate Y0 H0 Ys rhard];

    % Set strain history parameters
    edot          = 0.01; % strain rate
    emax          = 1.0;  % max strain
    num_reversals = 0;    % number of load reversals
    const2        = [edot emax num_reversals];

    % Call function to compute stress-strain history
    [e, sigma, t] = main_rd(const1, const2);

    % Plot the stress-strain curve
    figure;
    h=plot(e,sigma/1.E6,'k--');
    set(h,'LineWidth',2);
    switch pt
        case 1
            axis([0,1,0,600]);
        case 2
            axis([0,1,0,100]);
    end
    xlabel('$\epsilon$','Interpreter','latex');
    ylabel('$\sigma$, MPa','Interpreter','latex');
    set(gca,'FontSize',16)
    set(gca,'XMinorTick','On')
    set(gca,'YMinorTick','On')

    hold on

    % Set second strain rate
    edot          = 0.1; % strain rate
    emax          = 1.0; % max strain
    num_reversals = 0;   % number of load reversals
    const2        = [edot emax num_reversals];

    % Call function to compute stress-strain history
    [e, sigma, t] = main_rd(const1, const2);

    % Add to plot of stress-strain curve
    h=plot(e,sigma/1.E6,'k-.');
    set(h,'LineWidth',2);

    % Set third strain rate
    edot          = 1.0; % strain rate
    emax          = 1.0; % max strain
    num_reversals = 0;   % number of reversals
    const2        = [edot emax num_reversals];

    % Call function to compute stress-strain history
    [e, sigma, t] = main_rd(const1, const2);

    % Add to plot the stress-strain curve
    h=plot(e,sigma/1.E6,'k');
    set(h,'LineWidth',2);
    legend_handle = legend('$\dot\epsilon=0.01$/s',...
        '$\dot\epsilon=0.1$/s','$\dot\epsilon=1.0$/s');
    set(legend_handle,'FontSize',16,'Interpreter',...
        'latex','Location','Best')
    hold off
end

