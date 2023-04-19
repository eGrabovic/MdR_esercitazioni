function plot_joint_trajectories(qsol, options)
%
%
%
arguments
    qsol
    options.fontSize = 20;
    options.linewidth = 1.0;
    options.lineColor = 'k';
end

figure('color', 'w');

ax = subplot(3, 3, 1);
plot(ax, qsol(end,:), qsol(1,:), 'color', options.lineColor, 'linewidth', options.linewidth)
xlabel('t'); ylabel('$q_1$'); grid on
set(ax, 'Fontsize', options.fontSize);


ax = subplot(3, 3, 2);
plot(ax, qsol(end,:), qsol(1,:), 'color', options.lineColor, 'linewidth', options.linewidth)
grid on
set(ax, 'Fontsize', options.fontSize);
xlabel('t')
ylabel('$q_2$')

ax = subplot(3, 3, 3);
plot(ax, qsol(end,:), qsol(3,:), 'color', options.lineColor, 'linewidth', options.linewidth)
grid on
set(ax, 'Fontsize', options.fontSize);
xlabel('t')
ylabel('$q_3$')

ax = subplot(3, 3, 4);
plot(ax, qsol(end,:), qsol(4,:), 'color', options.lineColor, 'linewidth', options.linewidth)
grid on
set(ax, 'Fontsize', options.fontSize);
xlabel('t')
ylabel('$q_4$')

ax = subplot(3, 3, 5);
plot(ax, qsol(end,:), qsol(5,:), 'color', options.lineColor, 'linewidth', options.linewidth)
grid on
set(ax, 'Fontsize', options.fontSize);
xlabel('t')
ylabel('$q_5$')

ax = subplot(3, 3, 6);
plot(ax, qsol(end,:), qsol(6,:), 'color', options.lineColor, 'linewidth', options.linewidth)
grid on
set(ax, 'Fontsize', options.fontSize);
xlabel('t')
ylabel('$q_6$')

ax = subplot(3, 3, 7);
plot(ax, qsol(end,:), qsol(7,:), 'color', options.lineColor, 'linewidth', options.linewidth)
grid on
set(ax, 'Fontsize', options.fontSize);
xlabel('t')
ylabel('$q_7$')

ax = subplot(3, 3, 8);
plot(ax, qsol(end,:), qsol(8,:), 'color', options.lineColor, 'linewidth', options.linewidth)
grid on
set(ax, 'Fontsize', options.fontSize);
xlabel('t')
ylabel('$q_8$')

ax = subplot(3, 3, 9);
plot(ax, qsol(end,:), qsol(9,:), 'color', options.lineColor, 'linewidth', options.linewidth)
grid on
set(ax, 'Fontsize', options.fontSize);
xlabel('t')
ylabel('$q_9$')