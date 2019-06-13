function [lines] = ycPlot( ...,
    x_data, y_data, x_label, y_label, title_, leg_list, extra ...,
)
% 20190302
% Chao Yin
% last change 5/7

%% buttons
set(groot,'defaulttextinterpreter','latex','defaulttextFontWeight','bold');  
set(groot, 'defaultAxesTickLabelInterpreter','latex');  
set(groot, 'defaultColorbarTickLabelInterpreter','latex');  
set(groot, 'defaultLegendInterpreter','latex');
set(groot,'defaulttextFontSize',17)
set(groot,'defaultAxesFontSize',17)

%% 

%plot_number = size(y_data, 2);
line_type = '-';
colors = 0;
use_marker = false;
if nargin == 7
for extra_index = 1: length(extra)
    extra_item = extra{extra_index};
    if strcmp( extra_item{1}, 'Marker')
        use_marker = true;
    end
    if strcmp( extra_item{1}, 'line_type')
        line_type = extra_item{2};
    end
    if strcmp( extra_item{1}, 'color')
        colors = extra_item{2};
    end
end
end % if (nargin==7)

%% plot with error bar
err_plot = false;
if nargin == 7
for extra_index = 1: length(extra)
    extra_item = extra{extra_index};
    if strcmp( extra_item{1}, 'err')
        err_plot = true;
        hold on
        for curve_index = 1: size(y_data, 1)
            errorbar(x_data, y_data(curve_index, :), extra_item{2}(curve_index, :) ,'Linewidth',2 );
        end
    end
end % for extra_index = 1: length(extra)
end % if (nargin==7)

lines = [];
if ~err_plot & use_marker
    lines = plot( x_data, y_data, line_type, 'Linewidth',2  ,'Marker', '+', 'MarkerSize', 8);
elseif ~err_plot
    if colors == 0
        lines = plot( x_data, y_data, line_type, 'Linewidth', 1);
    else
        lines = plot( x_data, y_data, line_type, 'Linewidth', 1, 'color', colors);
    end
end
%loglog( x_data, y_data, 'Linewidth',2 ,'Marker', 'o', 'MarkerSize',8 );
%set( gca, 'yscale', 'log')

%% other lines

hold on
if nargin == 7
for extra_index = 1: length(extra)
    extra_item = extra{extra_index};
    if strcmp( extra_item{1}, 'x_par')
        y_lim=get(gca,'Ylim');
        period_line = plot( extra_item{2}*[1 1], y_lim,'m--','Linewidth',2 );
    end
    if strcmp( extra_item{1}, 'y_par')
        x_lim=get(gca,'Xlim');
        period_line = plot(x_lim, extra_item{2}*[1 1], 'm--','Linewidth',2 );
    end
end

end % if (nargin==7)

%% legend, label, title

leg = legend(leg_list);
%set(leg,'box','off' );% , 'Orientation','horizontal');
%title(leg, 'K')


xlabel(x_label)
ylabel(y_label)
title(title_)

end

