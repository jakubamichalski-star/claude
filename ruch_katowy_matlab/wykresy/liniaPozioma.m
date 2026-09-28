function liniaPozioma(ax, y, kolor, styl)
%LINIAPOZIOMA Pozioma linia odniesienia na całej szerokości osi (bez legendy).
if nargin < 3, kolor = [0.55 0.23 0.29]; end
if nargin < 4, styl = '--'; end
xl = xlim(ax);
plot(ax, xl, [y y], styl, 'Color', kolor, 'LineWidth', 0.8, 'HandleVisibility', 'off');
xlim(ax, xl);
end
