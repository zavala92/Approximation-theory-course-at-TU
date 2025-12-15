function padetables_demo
% PADETABLES_DEMO  Interactive Padé table demo.
% Author: Gentian Zavalani 14/12/2025

clc

% ---------- ask user for function and degrees ----------
fstr = input('f(z)? ','s');     % e.g. exp(z), cos(z), tan(z.^4)
f    = @(z) eval(fstr);

nmax = input('nmax? ');

% ---------- options ----------
printDegrees = false;    % console printing
drawInCells  = false;    % draw (mu,nu) labels on each colored square

% ---------- build colour matrix + store labels ----------
d      = zeros(nmax+1, nmax+1);               
labels = strings(nmax+1, nmax+1);

rand('state',7);                           % fixed seed (recommended)
h = tan(2*rand(nmax+1) - 1);                 % sized to nmax
h(min(20,nmax+1), 1) = 1;                  % keep your special entry safely

for n = 0:nmax
    for m = 0:nmax
        [~,~,~,mu,nu] = padeapprox(f,m,n);

        % store label for this cell
        labels(n+1,m+1) = sprintf('(%d,%d)', mu, nu);

        % Optional console output
        if printDegrees
            fprintf(' %s', labels(n+1,m+1));
        end

        % --- robust guard: mu,nu must be valid indices ---
        if ~isscalar(mu) || ~isscalar(nu) || ...
                ~isfinite(mu) || ~isfinite(nu) || ...
                mu < 0 || nu < 0 || ...
                mu ~= floor(mu) || nu ~= floor(nu) || ...
                (mu+1) > size(h,1) || (nu+1) > size(h,2)

            d(n+1,m+1) = NaN;
            labels(n+1,m+1) = "(fail)";
        else
            d(n+1,m+1) = h(mu+1,nu+1);
        end
    end
    if printDegrees
        fprintf('\n');
    end
end

% ---------- plot Padé table ----------
fig = figure(1); clf(fig);
set(fig,'Color','w','Renderer','painters');  % vector-safe
set(fig,'Units','inches','Position',[1 1 6.5 6.5]);

ax = axes(fig); %#ok<LAXES>

% Use imagesc so the full matrix is shown (pcolor drops last row/col)
imagesc(ax, 0:nmax, 0:nmax, d);

axis(ax,'ij');
axis(ax,'square');

set(ax, ...
    'XAxisLocation','top', ...
    'TickLength',[0 0], ...
    'XTick',0:nmax, ...
    'YTick',0:nmax, ...
    'FontSize',10);

colormap(ax, turbo);
set(ax,'Color',[1 1 1]);   % NaNs show as white

xlabel(ax,'m');
ylabel(ax,'n');
% title(ax,fstr,'Interpreter','tex');

% ---- optional: draw grid lines like the paper figure ----
hold(ax,'on');
for k = -0.5:1:(nmax+0.5)
    xline(ax,k,'k-','LineWidth',0.25);
    yline(ax,k,'k-','LineWidth',0.25);
end

% ---- draw (mu,nu) text labels centered in each cell ----
if drawInCells
    for n = 0:nmax
        for m = 0:nmax
            txt = labels(n+1,m+1);
            % use latex only if you want; here we keep it simple
            text(ax, m, n, txt, ...
                'Interpreter','none', ...
                'HorizontalAlignment','center', ...
                'VerticalAlignment','middle', ...
                'FontSize',10, ...
                'Color','k');
        end
    end
end

hold(ax,'off');
drawnow

% ---------- EXPORT PDF ----------
fname = regexprep(fstr,'[^a-zA-Z0-9]','');
if isempty(fname)
    fname = 'pade_table';
end
pdfname = ['padetable_' fname '.pdf'];

exportgraphics(fig,pdfname, ...
    'ContentType','vector', ...
    'BackgroundColor','white');

fprintf('Saved Padé table to %s\n', pdfname);

end
