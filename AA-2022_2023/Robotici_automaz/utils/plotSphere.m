function  plotSphere(r, C, options)
%
%
%
arguments
    r
    C
    options.parent = [];
    options.FaceColor = 'r';
    options.FaceAlpha = 0.5;
    options.Edgecolor = 'none';
    options.NFaces = 50;
end

[X,Y,Z] = sphere(options.NFaces);

X = X.*r+C(1);
Y = Y.*r+C(2);
Z = Z.*r+C(3);

if isempty(options.parent)
    figure; hold on;axis equal;
    options.parent = gca;
end

surface(X,Y,Z,'FaceColor',options.FaceColor, 'FaceAlpha', options.FaceAlpha, 'EdgeColor', options.Edgecolor, 'Parent', options.parent);

