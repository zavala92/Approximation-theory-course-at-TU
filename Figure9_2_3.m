%% AAA quadrature 
% Author: Gentian Zavalani 13.12.2025
clear; close all; clc;

% -----------------------------------------------------
%  Define the Bernstein ellipse Γ with foci ±1
% -----------------------------------------------------
rho = 1/sqrt(25) + sqrt(26/25);     % from the paper
c   = rho * exp(2i*pi*(1:200)/200); % sample on the ellipse
s   = (c + 1./c)/2;                 % map to ellipse in z-plane

% ep = 1/sqrt(25); 
% S = [-1i*ep+linspace(-1,1,100)'; 1+ep*exp(1i*pi*(-49:49)'/100)];
% s = [S; -S];

% -----------------------------------------------------
%  Compute C(s) = (2πi)^(-1) log((s+1)/(s-1))
% -----------------------------------------------------
C = log((s+1)./(s-1));

% -----------------------------------------------------
%  AAA rational approximation r_n ≈ 2πi*C
% -----------------------------------------------------
[r,pol,res] = aaa(C, s, 'degree', 20, 'sign', 1);

% --------- Style parameters (for polished plots) ---------
FS = 18;      % font size
LW = 2;       % line width
MS = 16;      % marker size

% -----------------------------------------------------
%  Plot AAA poles on the ellipse
% -----------------------------------------------------
figure('Color','w','Position',[100 100 1400 700]); % big, white figure

% Panel 1: AAA poles
subplot(2,2,1)
plot(real(s), imag(s), 'k-', 'LineWidth', LW); hold on;
plot(real(pol), imag(pol), 'ro', 'MarkerSize', MS, 'LineWidth', 1.5);
axis equal; box on;
title('AAA Poles','FontSize',FS+2)
xlabel('Re$(s)$','FontSize',FS)
ylabel('Im$(s)$','FontSize',FS)
set(gca,'FontSize',FS, 'LineWidth', 1, 'TickDir','out')

% Panel 2: Phase portrait of r(s)
subplot(2,2,3)
[xg,yg] = meshgrid(linspace(-1.2,1.2,500), linspace(-1,1,500));
zg = xg + 1i*yg;
rg = feval(r,zg);

% Phase portrait with a professional colormap
imagesc(linspace(-1.2,1.2,500), linspace(-1,1,500), angle(rg));
set(gca,'YDir','normal')
colormap(parula)   % More professional colormap
colorbar;
hold on;
plot(real(s), imag(s), 'k-', 'LineWidth', LW)
axis equal tight;
title('Phase Portrait of $r_{20}(z)$','FontSize',FS+2)
xlabel('Re$(z)$','FontSize',FS)
ylabel('Im$(z)$','FontSize',FS)
set(gca,'FontSize',FS)

% Panel 3: Quadrature error vs degree n
subplot(1,2,2)
f  = @(z) 1./(1+25*z.^2);        % test integrand
I_exact = quadgk(f,-1,1);        % high accuracy reference

nn = 2:2:80;                     % degrees
errAAA   = zeros(size(nn));
errGauss = zeros(size(nn));

for j = 1:length(nn)
    n = nn(j);

    % --- AAA quadrature ---
    [rj, polj, resj] = aaa(C, s, 'degree', n, 'sign', 1);
    I_AAA = sum(resj .* f(polj));
    errAAA(j) = abs(I_AAA - I_exact);

    % --- Gauss-Legendre quadrature ---
    [x,w] = legpts(n);           % Chebfun Gauss-Legendre points
    I_G = sum(w' .* f(x));
    errGauss(j) = abs(I_G - I_exact);
end

loglog(nn, errGauss, 'ko', 'MarkerSize', MS); hold on;
loglog(nn, errAAA,  'ro','MarkerSize',MS/2,'LineWidth',LW);
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
set(gca, 'XScale', 'linear', 'YScale', 'log')
xlim([0 80]); % x-axis from 0 to 80
exportgraphics(gcf,'aaa_quadrature.pdf', 'ContentType', 'vector', 'BackgroundColor', 'white');
