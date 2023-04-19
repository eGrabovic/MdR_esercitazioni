%% common scene options

% global frame
plotFrame(eye(4), 'scale', 0.5, 'label', '0', 'parent', gca);

% scene plot options
view(90,0)
xlabel('x')
ylabel('y')
zlabel('z')

set(gca, 'zlim', [-1.0 1.0]);
set(gca, 'xlim', [-1.0 1.0]);
set(gca, 'ylim', [-1.0 1.0]);
set(gca, 'box', 'on')