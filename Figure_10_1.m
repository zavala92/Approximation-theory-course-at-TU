function playpade
% PLAYPADE   Padé vs robust Padé on the unit disk (Chebfun, disk sampling)
%Author Gentian Zavalani 14/12/2025

clc

% ---------- ask user for function and degrees ----------
fstr = input('f(z)? ','s');      % e.g. 'log(1.2-z)' or 'tan(z.^4)'
f    = @(z) eval(fstr);

m = input('m? ');                % numerator degree
n = input('n? ');                % denominator degree

% ---------- parameters ----------
r  = 1;                          % sampling circle radius (< nearest sing.)
N  = 8*(m+n) + 1;                % oversampling on |z|=r

% error grid in |z| <= 0.5, odd multiples of 0.01
xs = -0.49:0.02:0.49;
ys = xs;
[X,Y] = meshgrid(xs,ys);
mask   = (X.^2 + Y.^2) <= 0.5^2;
zgrid  = X(mask) + 1i*Y(mask);
fgrid  = f(zgrid);

% circle to draw
th  = linspace(0,2*pi,600);
uc  = r*exp(1i*th);

% Make the circle look like your 2nd screenshot:
axlims = 1.4*[-1 1 -1 1];

% ===================== FIGURE STYLE (like playdisk) ======================
fig = figure(1); clf(fig);
set(fig,'Color','w','Name',fstr,'NumberTitle','off');
set(fig,'Units','inches','Position',[0.8 0.8 9.0 4.2]);
set(fig,'Renderer','painters'); % vector PDF

set(fig,'DefaultTextInterpreter','latex');
set(fig,'DefaultAxesTickLabelInterpreter','latex');
set(fig,'DefaultLegendInterpreter','latex');

tlo = tiledlayout(fig,1,2,'Padding','compact','TileSpacing','compact');

titleFont = 12;
textFont  = 10;
circleLW  = 2.2;
boxLW     = 1.2;
poleMS    = 8;

%% ================= 1. Padé, tol = 0 ==========================
tol0 = 0;
tic
[r0,a0,b0,mu0,nu0,poles0,res0] = padeapprox(f,m,n,tol0,r,N);
t0   = toc;

format long g
disp('Zeros:'), disp(roots(a0(end:-1:1)))
disp('Poles:'), disp(roots(b0(end:-1:1)))

rgrid0 = r0(zgrid);
err0   = norm(fgrid - rgrid0, inf);

ax = nexttile(tlo,1); hold(ax,'on');

% thick circle (NO fill, NO dotted)
plot(ax,real(uc),imag(uc),'k-','LineWidth',circleLW);

% poles colored by |residue|
absRes0 = abs(res0);
for k = 1:numel(poles0)
    c = residue_color(absRes0(k));
    plot(ax,real(poles0(k)),imag(poles0(k)), ...
        'o','MarkerFaceColor',c,'MarkerEdgeColor','k','MarkerSize',poleMS);
end

axis(ax,'equal');

% small padding so box not clipped
pad = 0.02;
xl = axlims(1:2); yl = axlims(3:4);
dx = diff(xl);    dy = diff(yl);
axis(ax,[xl(1)-pad*dx, xl(2)+pad*dx, yl(1)-pad*dy, yl(2)+pad*dy]);

set(ax,'XTick',[],'YTick',[], ...
    'Box','on','LineWidth',boxLW,'FontSize',textFont);

title(ax,'Pad\''e (tol = 0)', 'FontSize',titleFont,'FontWeight','normal');


text(ax,0.97,0.96,sprintf('$(%d,%d)$',m,n), ...
    'Units','normalized','HorizontalAlignment','right', ...
    'VerticalAlignment','top','FontSize',textFont,'Clipping','off');

text(ax,0.03,0.06,sprintf('$(%d,%d)$\n$%.3f\\,\\mathrm{s}$',mu0,nu0,t0), ...
    'Units','normalized','HorizontalAlignment','left', ...
    'VerticalAlignment','bottom','FontSize',textFont,'Clipping','off');

text(ax,0.97,0.06,sprintf('$\\mathrm{Err}=%.2e$',err0), ...
    'Units','normalized','HorizontalAlignment','right', ...
    'VerticalAlignment','bottom','FontSize',textFont,'Clipping','off');

%% ================= 2. Robust Padé ============================
tolR = 1e-14;
tic
[r1,~,~,mu1,nu1,poles1,res1] = padeapprox(f,m,n,tolR,r,N); %#ok<ASGLU>
t1   = toc;

rgrid1 = r1(zgrid);
err1   = norm(fgrid - rgrid1, inf);

ax = nexttile(tlo,2); hold(ax,'on');

plot(ax,real(uc),imag(uc),'k-','LineWidth',circleLW);

absRes1 = abs(res1);
for k = 1:numel(poles1)
    c = residue_color(absRes1(k));
    plot(ax,real(poles1(k)),imag(poles1(k)), ...
        'o','MarkerFaceColor',c,'MarkerEdgeColor','k','MarkerSize',poleMS);
end

axis(ax,'equal');
pad = 0.02;
xl = axlims(1:2); yl = axlims(3:4);
dx = diff(xl);    dy = diff(yl);
axis(ax,[xl(1)-pad*dx, xl(2)+pad*dx, yl(1)-pad*dy, yl(2)+pad*dy]);

set(ax,'XTick',[],'YTick',[], ...
    'Box','on','LineWidth',boxLW,'FontSize',textFont);

title(ax,'Robust Pad\''e', 'FontSize',titleFont,'FontWeight','normal');


text(ax,0.97,0.96,sprintf('$(%d,%d)$',m,n), ...
    'Units','normalized','HorizontalAlignment','right', ...
    'VerticalAlignment','top','FontSize',textFont,'Clipping','off');

text(ax,0.03,0.06,sprintf('$(%d,%d)$\n$%.3f\\,\\mathrm{s}$',mu1,nu1,t1), ...
    'Units','normalized','HorizontalAlignment','left', ...
    'VerticalAlignment','bottom','FontSize',textFont,'Clipping','off');

text(ax,0.97,0.06,sprintf('$\\mathrm{Err}=%.2e$',err1), ...
    'Units','normalized','HorizontalAlignment','right', ...
    'VerticalAlignment','bottom','FontSize',textFont,'Clipping','off');

drawnow;

% ===================== EXPORT (vector PDF) ==============================
exportgraphics(fig,'forisart.pdf','ContentType','vector','BackgroundColor','white');

end

% ------ |residue| → color (same scheme as playdisk) ----------
function c = residue_color(r)
% Color by |residue| exactly as in Trefethen–Gutknecht
if r >= 1e-3
    c = [0 0 1];          % blue: true poles
elseif r >= 1e-6
    c = [0.3 0.6 1.0];    % light blue
elseif r >= 1e-9
    c = [0 0.6 0];        % green
elseif r >= 1e-12
    c = [0.6 1.0 0.6];    % light green
elseif r >= 1e-14
    c = [1.0 0.6 0.8];    % pink: Froissart doublets
else
    c = [1 0 0];          % red: numerical junk
end
end

