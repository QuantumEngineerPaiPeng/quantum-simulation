function ppStyle(figure_FontSize,figure_LineWidth,figure_MarkerSize)

set(get(gca,'XLabel'),'FontSize',figure_FontSize,'Vertical','top');
set(get(gca,'YLabel'),'FontSize',figure_FontSize,'Vertical','middle');
set(findobj('FontSize',10),'FontSize',figure_FontSize);
if nargin>1
    set(findobj(get(gca,'Children'),'LineWidth',0.5),'LineWidth',figure_LineWidth);
    if nargin>2
        set(findobj(get(gca,'Children'),'MarkerSize',6),'MarkerSize',figure_MarkerSize);
    end
end
end