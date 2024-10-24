function [transforms, ax] = init_plot_DH(Tj, blist, Jtype_list, options)
%
%
%

arguments
    Tj
    blist
    Jtype_list
    options.joint_len = 0.5;
    options.joint_r = 0.2;
    options.T0 = eye(4);
    options.TE = eye(4);
    options.parent = [];
    options.label = [];
end

    function [ptX, ptY, ptZ] = createCylinder(r, lentop, lenbot, axis, axialOff)
        
        theta = linspace(0, 2*pi, 15);
        x = r.*cos(theta);
        y = r.*sin(theta);
        ztop = lentop.*ones(1,15);
        zbot = -lenbot.*ones(1,15);
        pttop = [x;y;ztop];
        ptbot = [x;y;zbot];
        theta = acos(axis(3));
        phi = atan2(axis(1), axis(2));
        R = rotZ(-phi)*rotX(-theta);
        pttop = R*pttop;
        ptbot = R*ptbot;
        
        ptX = [pttop(1,:); ptbot(1,:)];
        ptY = [pttop(2,:); ptbot(2,:)];
        ptZ = [pttop(3,:); ptbot(3,:)]-axialOff;
        
    end

    function [points, faces] = createParallelepiped(edge, len, axis, axialOff)
        
        a = edge;
        if length(len) == 2
            h1 = len(1);
            h2 = len(2);
        else
            h1 = 0;
            h2 = len;
        end
        P1 = [a/2;-a/2;-h1]; P2 = [a/2;a/2;-h1]; P3 = [a/2;a/2;h2]; P4 = [a/2;-a/2;h2];
        P5 = [-a/2;a/2;-h1]; P6 = [-a/2;a/2;h2]; P7 = [-a/2;-a/2;-h1]; P8 = [-a/2;-a/2;h2];
        points = [P1,P2,P3,P4,P5,P6,P7,P8] - [0; 0; axialOff];
        theta = acos(axis(3));
        phi = atan2(axis(1), axis(2));
        R = rotZ(-phi)*rotX(-theta);
        points = R*points;
        faces = [1,2,3,4;2,5,6,3; 3,6,8,4; 8,6,5,7; 1,4,8,7; 1,7,5,2];
        
    end

if isempty(options.parent)
    figure('color', 'w', 'units', 'normalized', 'outerposition', [0 0 1 1]);
    hold on;
    axis equal;
    ax = gca;  % get current axes 
else
    ax = options.parent;
end

nj = length(Tj);

transforms = cell(1, nj);
j_len = options.joint_len./2;
j_R = options.joint_r;
T = options.T0;
transform0 = hgtransform('parent', ax);
transform0.Matrix = options.T0;
parent = transform0;

% plot base frame
plotFrame(eye(4), 'parent', parent, 'label', sprintf('base%s', options.label), 'scale', 0.15);

for j = 1:nj
    
    transforms{j} = hgtransform('parent', ax);
    parentPrev = parent;
    parent = transforms{j};
    T = T*Tj{j};
    transforms{j}.Matrix = T;
    P = rigidInverse(Tj{j}); % O_i-1O_i vector
    P = P(1:3, 4);
    
    
    % geometric primitives
    if strcmpi(Jtype_list(j) , 'R')
        [jointX, jointY, jointZ] = createCylinder(j_R, j_len, j_len, [0;0;1], blist(j)); % joint primitive
        % plot joint
        surf(jointX, jointY, jointZ, 'facecolor', 'r', 'edgecolor', 'none', 'parent', parentPrev);
        fill3(jointX(1,:), jointY(1,:),  jointZ(1,:), 'r', 'parent', parentPrev);
        fill3(jointX(2,:), jointY(2,:),  jointZ(2,:), 'r', 'parent', parentPrev);
        
        % build link primitive
        [linkX, linkY, linkZ] = createCylinder(j_R/1.75, norm(P), 0, (P)./norm(P),0); % link primitive
        
    else
        [ptsP, faces] = createParallelepiped(j_R*2, [j_len*1.5 j_len*1.5], [0;0;1], blist(j));
        [ptsP_mobile, faces_mobile] = createParallelepiped(j_R*2*0.5, [j_len*1.5*3 j_len*1.5*0], [0;0;1], blist(j));
        patch('faces', faces, 'vertices', ptsP.',  'facecolor', 'green', 'Parent', parentPrev);
        patch('faces', faces_mobile, 'vertices', ptsP_mobile.',  'facecolor', 'green', 'Parent', parent);
        
         % build link primitive
        [linkX, linkY, linkZ] = createCylinder(j_R/1.75, j_len*10, 0, (P)./norm(P),0); % link primitive
        
        
    end
    % plot link
    surf(linkX, linkY, linkZ, 'facecolor', 'b', 'edgecolor', 'none', 'parent', parent);
    fill3(linkX(1,:), linkY(1,:),  linkZ(1,:), 'b', 'parent', parent);
    fill3(linkX(2,:), linkY(2,:),  linkZ(2,:), 'b', 'parent', parent);
end
transforms{end+1} = hgtransform('parent', parent);
transforms{end}.Matrix = options.TE;

% plot end effector
[EEX, EEY, EEZ] = createCylinder(j_R/2, j_len, j_len, [0;0;1], 0);

surf(EEX, EEY, EEZ, 'facecolor', '#4DBEEE', 'edgecolor', 'none', 'parent', transforms{end});
fill3(EEX(1,:), EEY(1,:), EEZ(1,:), EEZ(1,:).*0, 'facecolor', '#4DBEEE', 'parent', transforms{end});
fill3(EEX(2,:), EEY(2,:), EEZ(2,:), EEZ(1,:).*0, 'facecolor', '#4DBEEE', 'parent', transforms{end});

plotFrame(eye(4), 'parent', transforms{end}, 'scale', 0.1, 'label', sprintf('EE%s', options.label))
end