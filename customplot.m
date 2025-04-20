function customplot(plt,xtext,ytext)

ax = gca;

% Lines and colors
plt.LineWidth = 2;
% plt.LineColor = [0.07 0.62 1 0.8];

% Extend the limits a little
y = get(ax,'YLim');
y(1) = y(1)-0.1*y(1);
y(2) = y(2)+0.1*y(2);
ax.YLim = y;

% Set fonts
ax.FontName = 'Times New Roman';
ax.FontSize = 16;

% Set ticks and labels
ax.XTickLabel = strrep(xticklabels,'-',char(8722));
ax.YTickLabel = strrep(yticklabels,'-',char(8722));
ax.TickDir = 'out';
ax.LineWidth = 1.2;
ax.XLabel.String = xtext;
ax.YLabel.String = ytext;

% Grids and boxes
grid on
box off
