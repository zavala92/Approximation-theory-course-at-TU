function playdisk
% PLAYDISK  Rational interpolation on the unit disk using RATDISK.
% 2x2 figure: Interpolation / Least-squares, unstabilized vs stabilized (SVD).
% Fig-style annotations:
%   top-right: (m,n,N)
%   bottom-left: (mu,nu) and time
%   bottom-right: Err
% Also shows sampling nodes (N+1 roots of unity) and the unit circle.
% Exports vector PDF: forisart.pdf
% Author: Gentian Zavalani 13.12.2025

clc

% ---- Ask user for function and degrees ----
fstr = input('f(z)? ','s');      % e.g. 'exp(z)./(1.1-z.^2)' or 'log(2+z.^4)./(1-16*z.^4)'
f    = @(z) eval(fstr);          % uses variable z in fstr

m = input('m? ');                % numerator degree
n = input('n? ');                % denominator degree

% ---------------- parameters ----------------
N1  = m + n;                     % interpolation: N1+1 sample points
N2  = 4*(m+n) + 1;               % LS: oversampled
tol = 1e-14;                     % SVD stabilization tolerance

% ---- grid in the unit disk ----
xs = -0.99:0.02:0.99;
ys = xs;
[X,Y] = meshgrid(xs,ys);
mask   = (X.^2 + Y.^2) <= 1;
zgrid  = X(mask) + 1i*Y(mask);
fgrid  = f(zgrid);

% ---- unit circle for boundary ----
theta  = linspace(0,2*pi,400);
ucirc  = exp(1i*theta);
axlims = 1.4*[-1 1 -1 1];

% ===================== COMPUTE (4 fits) ================================
[tI,  muI,  nuI,  errI,  polesI,  absResI,  NusedI ] = one_fit(f,m,n,N1,0,  zgrid,fgrid);
[tLS, muLS, nuLS, errLS, polesLS, absResLS, NusedLS] = one_fit(f,m,n,N2,0,  zgrid,fgrid);
[tIS, muIS, nuIS, errIS, polesIS, absResIS, NusedIS] = one_fit(f,m,n,N1,tol,zgrid,fgrid);
[tLSS,muLSS,nuLSS,errLSS,polesLSS,absResLSS,NusedLSS]= one_fit(f,m,n,N2,tol,zgrid,fgrid);

% ===================== PLOTTING========================
fig = figure(1); clf(fig);
set(fig,'Color','w','Name',fstr,'NumberTitle','off');
set(fig,'Units','inches','Position',[0.8 0.8 9.0 7.0]);
set(fig,'Renderer','painters'); % vector PDF

% LaTeX everywhere (MATLAB subset)
set(fig,'DefaultTextInterpreter','latex');
set(fig,'DefaultAxesTickLabelInterpreter','latex');
set(fig,'DefaultLegendInterpreter','latex');

tlo = tiledlayout(fig,2,2,'Padding','compact','TileSpacing','compact');

titleFont = 12;
textFont  = 10;
circleLW  = 1.6;
nodeMS    = 5;
boxLW     = 1.0;
poleMS    = 7;

% --- panel 1: interpolation (tol=0)
ax = nexttile(tlo,1);
plot_panel(ax,ucirc,NusedI,polesI,absResI,axlims, ...
    'Interpolation (tol = 0)', m,n,NusedI,muI,nuI,tI,errI, ...
    circleLW,nodeMS,boxLW,poleMS,titleFont,textFont);

% --- panel 2: least-squares (tol=0)
ax = nexttile(tlo,2);
plot_panel(ax,ucirc,NusedLS,polesLS,absResLS,axlims, ...
    'Least-squares (tol = 0)', m,n,NusedLS,muLS,nuLS,tLS,errLS, ...
    circleLW,nodeMS,boxLW,poleMS,titleFont,textFont);

% --- panel 3: interpolation (stabilized)
ax = nexttile(tlo,3);
plot_panel(ax,ucirc,NusedIS,polesIS,absResIS,axlims, ...
    'Interpolation (stabilized)', m,n,NusedIS,muIS,nuIS,tIS,errIS, ...
    circleLW,nodeMS,boxLW,poleMS,titleFont,textFont);

% --- panel 4: least-squares (stabilized)
ax = nexttile(tlo,4);
plot_panel(ax,ucirc,NusedLSS,polesLSS,absResLSS,axlims, ...
    'Least-squares (stabilized)', m,n,NusedLSS,muLSS,nuLSS,tLSS,errLSS, ...
    circleLW,nodeMS,boxLW,poleMS,titleFont,textFont);

drawnow;

% ===================== EXPORT (vector PDF) ==============================
exportgraphics(fig,'forisart.pdf','ContentType','vector','BackgroundColor','white');

end

%% ===================== numerical helper ============================
function [t_elapsed,mu,nu,Err,poles,absRes,Nused] = one_fit(f,m,n,N,tol,zgrid,fgrid)
tstart = tic;
[r,~,~,mu,nu,poles,residues] = ratdisk(f, m, n, N, tol); %#ok<ASGLU>
t_elapsed = toc(tstart);

rgrid  = r(zgrid);
Err    = norm(fgrid - rgrid, inf);
absRes = abs(residues);

% N is the sampling parameter (N+1 roots of unity)
Nused = N;
end

%% ===================== plotting helper (Fig-style) ======================
function plot_panel(ax,ucirc,N,poles,absRes,axlims, ...
                    mode,m,n,Ndisp,mu,nu,t_elapsed,Err, ...
                    circleLW,nodeMS,boxLW,poleMS,titleFont,textFont)

cla(ax); hold(ax,'on');

% unit circle
plot(ax, real(ucirc), imag(ucirc), 'k-', 'LineWidth', circleLW);

% sampling nodes (N+1 roots of unity)
zj = exp(2i*pi*(0:N)/(N+1));
plot(ax, real(zj), imag(zj), 'k.', 'MarkerSize', nodeMS);

% poles colored by |residue|
for k = 1:numel(poles)
    c = residue_color(absRes(k));
    plot(ax, real(poles(k)), imag(poles(k)), 'o', ...
        'MarkerFaceColor', c, ...
        'MarkerEdgeColor', 'k', ...
        'MarkerSize', poleMS);
end

axis(ax,'equal');

% ---- FIX: small padding so the box is never clipped by tiledlayout ----
pad = 0.02; % 2% margin
xl = axlims(1:2); yl = axlims(3:4);
dx = diff(xl);    dy = diff(yl);
axis(ax,[xl(1)-pad*dx, xl(2)+pad*dx, yl(1)-pad*dy, yl(2)+pad*dy]);

set(ax,'XTick',[],'YTick',[], ...
    'Box','on', ...
    'LineWidth',boxLW, ...
    'FontSize',textFont);

title(ax, ['\textbf{' mode '}'], 'FontSize', titleFont, 'FontWeight','normal');

% (m,n,N) top-right
text(ax,0.97,0.96,sprintf('$(%d,%d,%d)$',m,n,Ndisp), ...
    'Units','normalized','HorizontalAlignment','right', ...
    'VerticalAlignment','top','FontSize',textFont,'Clipping','off');

% (mu,nu) and time bottom-left
text(ax,0.03,0.06,sprintf('$(%d,%d)$\n$%.3f\\,\\mathrm{s}$',mu,nu,t_elapsed), ...
    'Units','normalized','HorizontalAlignment','left', ...
    'VerticalAlignment','bottom','FontSize',textFont,'Clipping','off');

% Err bottom-right
text(ax,0.97,0.06,sprintf('$\\mathrm{Err}=%.2e$',Err), ...
    'Units','normalized','HorizontalAlignment','right', ...
    'VerticalAlignment','bottom','FontSize',textFont,'Clipping','off');

end

%% ===================== residue color mapping ============================
function c = residue_color(r)
if r >= 1e-3
    c = [0 0 1];               % blue
elseif r >= 1e-6
    c = [0.4 0.7 1.0];         % light blue
elseif r >= 1e-9
    c = [0 0.6 0];             % green
elseif r >= 1e-12
    c = [0.6 1.0 0.6];         % light green
elseif r >= 1e-14
    c = [1.0 0.6 0.8];         % pink
else
    c = [1 0 0];               % red
end
end
