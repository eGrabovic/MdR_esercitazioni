function plot_orientation_errors(e_num, options)
%
%
%
arguments
    e_num
    options.fontSize = 20;
    options.linewidth = 1.0;
    options.lineColor = 'k';
end

figure('color', 'w');

ax = subplot(2, 3, 1);
plot(ax, e_num(1,:), 'color', options.lineColor, 'linewidth', options.linewidth)
xlabel('t'); ylabel('$e_x$'); grid on
set(ax, 'Fontsize', options.fontSize);


ax = subplot(2, 3, 2);
plot(ax, e_num(1,:), 'color', options.lineColor, 'linewidth', options.linewidth)
grid on
set(ax, 'Fontsize', options.fontSize);
xlabel('t')
ylabel('$e_y$')

ax = subplot(2, 3, 3);
plot(ax, e_num(3,:), 'color', options.lineColor, 'linewidth', options.linewidth)
grid on
set(ax, 'Fontsize', options.fontSize);
xlabel('t')
ylabel('$e_z$')

ax = subplot(2, 3, 4);
plot(ax, e_num(4,:), 'color', options.lineColor, 'linewidth', options.linewidth)
grid on
set(ax, 'Fontsize', options.fontSize);
xlabel('t')
ylabel('$e_{o1}$')

ax = subplot(2, 3, 5);
plot(ax, e_num(5,:), 'color', options.lineColor, 'linewidth', options.linewidth)
grid on
set(ax, 'Fontsize', options.fontSize);
xlabel('t')
ylabel('$e_{o2}$')

ax = subplot(2, 3, 6);
plot(ax, e_num(6,:), 'color', options.lineColor, 'linewidth', options.linewidth)
grid on
set(ax, 'Fontsize', options.fontSize);
xlabel('t')
ylabel('$e_{o3}$')