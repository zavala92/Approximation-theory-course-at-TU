 %% Rational Interpolation Demo for Lecture
  % Author: Gentian Zavalani 13.12.2025
% f(x) = exp(x) + 0.3 sin(3x) on [-1,1]

clear; close all;

% Setup LaTeX rendering
set(0,'defaultTextInterpreter','latex');
set(0,'defaultAxesTickLabelInterpreter','latex');
set(0,'defaultLegendInterpreter','latex');

% Function
% f = chebfun(@(x) cos(exp(x)), [-1 1]);
f = chebfun(@(x) exp(x) + 0.3*sin(3*x), [-1 1]);

colors = lines(6);    % MATLAB's distinct color set
for n = 1:7
    % ----- Compute rational (n,n) interpolant
    [p, q] = ratinterp(f, n, n);
    r = p ./ q;
    
    % ----- Chebyshev nodes (2n+1)
    k  = 0:(2*n);
    xk = cos(pi*k/(2*n));
    fk = f(xk);
    
    % ----- Plot
    figure(1); clf;
    set(gcf,'Color','w');
    hold on; box on; grid on;
    
    % Plot original function (black)
    plot(f, 'Color', [0 0 0], 'LineWidth', 2.2);
    
    % Plot rational interpolant (color varies with n)
    plot(r, '-', 'Color', colors(2,:), 'LineWidth', 2.5);
    
    % Chebyshev nodes
    plot(xk, fk, 'o', ...
         'MarkerSize', 8, ...
         'MarkerFaceColor', colors(2,:), ...
         'MarkerEdgeColor', 'k', ...
         'LineWidth', 1.2);
    
    % ----- Error
    err = norm(f - r, inf);
    
    % ----- Title (proper LaTeX)
    title(sprintf('$n = %d \\quad \\mathrm{error} = %7.2e$', n, err), ...
          'FontSize', 20, 'Interpreter', 'latex');
    
    xlabel('$x$', 'FontSize', 18);
    ylabel('$f(x)$, $r_n(x)$', 'FontSize', 18);
    
    axis([-1 1 0 3.2])
    set(gca,'FontSize',16,'LineWidth',1.3)
    
    % Legend (shows color for r_n)
    legend({'$f(x)$', sprintf('$r_{%d}(x)$',n)}, ...
           'Location','NorthWest', 'FontSize', 14);
    
    drawnow;
    disp('Press any key for next n...');
    pause;

end
