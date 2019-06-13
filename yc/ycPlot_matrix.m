function [outputArg1,outputArg2] = ycPlot_matrix( ...,
    plot_ma, x_label, y_label, title_, extra ...,
)
% 20190403
% Chao Yin
% last change 4/3

%% buttons
set(groot,'defaulttextinterpreter','latex','defaulttextFontWeight','bold');  
set(groot, 'defaultAxesTickLabelInterpreter','latex');  
set(groot, 'defaultColorbarTickLabelInterpreter','latex');  
set(groot, 'defaultLegendInterpreter','latex');
set(groot,'defaulttextFontSize',17)
set(groot,'defaultAxesFontSize',17)

%% pcolor
colorplot = pcolor( enlarge_pcolor( plot_ma ) );
set(colorplot , 'LineStyle','none');

%% color setting
%caxis([0 max(max(plot_ma)/2 )])
%colormap(jet)
%colorbar('Ticks',[0:0.002:0.006])
colorbar()

%{
x_n = size(plot_ma, 2);
xticks(linspace(1, x_n+1, 4) )
xticklabels({0,1,2,3})
%}

%% other lines
hold on
if nargin == 5
for extra_index = 1: length(extra)
    extra_item = extra{extra_index};
    if strcmp( extra_item{1}, 'x_par')
        y_lim=get(gca,'Ylim');
        period_line = plot( extra_item{2}*[1 1], y_lim,'r--','Linewidth', 2 );
    end
end

%% text
xlabel(x_label)
ylabel(y_label)
title(title_)

%% other lines
hold on
if nargin == 7
for extra_index = 1: length(extra)
    extra_item = extra{extra_index};
    if strcmp( extra_item{1}, 'x_par')
        y_lim=get(gca,'Ylim');
        period_line = plot( extra_item{2}*[1 1], y_lim,'m--','Linewidth',2 );
    end
end

end % if (nargin==7)

%% legend, label, title

if (nargin>=6)  
leg = legend(leg_list);
%set(leg,'box','off' );% , 'Orientation','horizontal');
%title(leg, 'K')

end % if (nargin==6)  

xlabel(x_label)
ylabel(y_label)
title(title_)
end

