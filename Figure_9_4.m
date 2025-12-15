%% ============================================================
% Application 2: Quadrature with Nearby Singularities
% f(z) = 1/(1 + 100 z^2)
% Author: Gentian Zavalani 13.12.2025
%% ============================================================

clear; close all; clc;

f = @(z) 1./(1 + 100*z.^2);

% ============================================================
% Step 1. Define Bernstein ellipse Γ with added slits
% ============================================================
rho   = 1.6;
theta = linspace(0,2*pi,400).';
c     = rho * exp(1i*theta);
s     = 0.5*(c + 1./c);                 % Bernstein ellipse (column)

% Add small vertical slits (for geometry illustration)
slit_len = 0.06;      % half-length of slit (choose small)
slit_n   = 80;

z1 =  1i/10;          % +0.1i singularity
z2 = -1i/10;          % -0.1i singularity

slit1 = z1 + 1i*linspace(-slit_len, slit_len, slit_n).';
slit2 = z2 + 1i*linspace(-slit_len, slit_len, slit_n).';

Gamma = [s; slit1; slit2];              % IMPORTANT: column vector for AAA

% Cauchy transform samples
C = log((Gamma+1)./(Gamma-1));

% ============================================================
% Step 2. AAA approximation (base degree 20)
% ============================================================
deg0 = 20;
[r, pol, res] = aaa(C, Gamma, 'degree', deg0, 'sign', 1);

% ============================================================
% Step 3. Exact integral and Gauss–Legendre comparison
% ============================================================
I_exact = quadgk(f, -1, 1);

nn = 2:2:80;
errAAA = zeros(size(nn));
errGauss = zeros(size(nn));

for j = 1:length(nn)
    n = nn(j);

    % AAA with degree n
    [rj, polj, resj] = aaa(C, Gamma, 'degree', n, 'sign', 1);
    I_AAAj = sum(resj .* f(polj));
    errAAA(j) = abs(I_AAAj - I_exact);

    % Gauss-Legendre
    [x, w] = legpts(n);
    Ig = sum(w'.*f(x));
    errGauss(j) = abs(Ig - I_exact);
end

% ============================================================
% Professional Plot Styling (same format as your example)
% ============================================================
FS = 18;     % font size
LW = 2;      % line width
MS = 16;     % marker size

figure('Color','w','Position',[100 100 1400 700]); % big, white figure

% ============================================================
% (1) AAA poles on ellipse with slits
% ============================================================
subplot(2,2,1)
plot(real(s), imag(s), 'k-', 'LineWidth', LW); hold on;
plot(real(slit1), imag(slit1), 'k-', 'LineWidth', LW);
plot(real(slit2), imag(slit2), 'k-', 'LineWidth', LW);

plot(real(pol), imag(pol), 'ro', 'MarkerSize', MS, 'LineWidth', 1.5);

axis equal; box on;
title('AAA Poles','FontSize',FS+2)
xlabel('Re$(s)$','FontSize',FS,'Interpreter','latex')
ylabel('Im$(s)$','FontSize',FS,'Interpreter','latex')
set(gca,'FontSize',FS, 'LineWidth', 1, 'TickDir','out')

% ============================================================
% (2) Phase portrait of AAA approximant
% ============================================================
subplot(2,2,3)

[xg,yg] = meshgrid(linspace(-1.5,1.5,600), linspace(-1,1,600));
zg = xg + 1i*yg;
rg = feval(r, zg);

imagesc(linspace(-1.5,1.5,600), linspace(-1,1,600), angle(rg));
set(gca,'YDir','normal')
colormap(parula)      % professional colormap (instead of hsv)

cb = colorbar;        % keep like your example (simple)
cb.TickLabelInterpreter = 'latex';

hold on;
plot(real(s), imag(s), 'k-', 'LineWidth', LW)
plot(real(slit1), imag(slit1), 'k-', 'LineWidth', LW)
plot(real(slit2), imag(slit2), 'k-', 'LineWidth', LW)

axis equal tight;
title(sprintf('Phase Portrait of $r_{%d}(z)$',deg0), ...
      'FontSize',FS+2,'Interpreter','latex')
xlabel('Re$(z)$','FontSize',FS,'Interpreter','latex')
ylabel('Im$(z)$','FontSize',FS,'Interpreter','latex')
set(gca,'FontSize',FS, 'LineWidth', 1, 'TickDir','out')

% ============================================================
% (3) Convergence plot: Gauss vs AAA
% ============================================================
subplot(1,2,2)

loglog(nn, errGauss, 'ko', 'MarkerSize', MS); hold on;
loglog(nn, errAAA,  'ro', 'MarkerSize', MS/2, 'LineWidth', LW);
grid on;

xlabel('Degree $n$','FontSize',FS,'Interpreter','latex')
ylabel('$|I_n - I|$','FontSize',FS,'Interpreter','latex')

legend({'Gauss-Legendre','AAA'}, ...
       'Location','southwest', ...
       'FontSize',FS, ...
       'Interpreter','latex')

title('\textbf{Convergence of quadrature}', ...
      'FontSize',FS+2, ...
      'Interpreter','latex')

set(gca,'FontSize',FS, 'LineWidth', 1, 'TickDir','out');
set(gca,'XScale','linear','YScale','log')
xlim([0 80])

% ============================================================
% Export (vector PDF for LaTeX)
% ============================================================
exportgraphics(gcf,'aaa_quadrature_nearby_singularities.pdf', ...
    'ContentType','vector','BackgroundColor','white');
