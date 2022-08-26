function [] = makeScaleBar(xlim,ylim,x_str_label,y_str_label)
% Make a scale bar, and take off plotting axes
x_al = xlim;
y_al = ylim;

plot([x_al(1); x_al(1)], [y_al(1); y_al(1)+(.25*range(y_al))], '-k',  [x_al(1);x_al(1)+(.25*range(x_al))], ...
      [y_al(1); y_al(1)], '-k', 'LineWidth', 2)
hold on
%axis([[10  34]    -60  140]) % move in scale bar
text((.25*range(x_al))/2,y_al(1)-.1*range(y_al), [num2str(range(x_al)*.25) ' ' x_str_label], 'HorizontalAlignment','center')
text((.25*range(x_al))/2,y_al(1)+.1*range(y_al), [num2str(range(x_al)*.25) ' ' y_str_label], 'HorizontalAlignment','left')
set(gca, 'Visible', 'off')
end