function  [ax, transform] = init_plot_body(G0)
%
% % G is the homogeneous transformation matrix
%

% to animate our body we employ the built-in matlab function "hgtransform"

%% Figure and axis
figure('color', 'w'); hold on; axis equal
title('Body simulation scene')
xlabel('x')
ylabel('y')
zlabel('z')
ax = gca; % "gca stands for GetCurrentAxes" if there are not any, Matlab will instantiate a figure with default axes
box(ax, 'on');
set(ax, 'Fontsize', 30);
grid(ax, 'on');
view(ax, 40, 35);
xlim(ax, [-5, 5]);
ylim(ax, [-5, 5]);
zlim(ax, [-4, 4]);

%% transform handle to update while animating
transform = hgtransform(ax);
transform.Matrix = G0;

% plot object graphics inside the transform 
[Faces, Vertices] = stlread('body.STL'); % read mesh file (note we overloaded the built-in matlab "stlread" function with another one inside "Functions" folder)
Vertices = Vertices+[-4,-0.75,-3+0.41]; % reposition the vertices to make the graphics center w.r.t. G
patch('Faces', Faces, 'Vertices', Vertices, 'edgecolor', 'none', 'faceColor', 'r', 'parent', transform);
light