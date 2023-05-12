function [ax, transforms] = plotRobot_localPOE(G0E, G0, X, q, Jtype, options)
% [ax, transforms] = plotRobot_localPOE(Glocal, G0, X, q)
% starts a graphics scene with a serial robot manipulator 
%
% INPUTS: -)Glocal: casadi function that defines relative homogeneous transformation
%                  between two links performed by each joint
%
%         -)G0: 1xN cell of 4x4 homogenous transofrmation matrices
%               representing the offset (relative) transformations between the
%               links (Glocal(zeros))
%
%         -) X: unitary distal twists
%
%         -) q: Nx1 of joint numerical values
%
% OUTPUTS:
%         -) ax: handle for the axis
%         -) transforms: handle for hgtransforms (matlab graphical tool to
%              manage homogeneous transformations between graphical objects)

arguments
    G0E
    G0
    X
    q
    Jtype
    options.jointradius = 0.03;
    options.jointlength = 0.08;
end

figure; hold on; axis equal; box on; grid on;
ax = gca;  % get current axes

G = eye(4);
transforms = cell(1, length(q));
Glocalnum = cell(1, length(q));
[Glocalnum{:}] = G0E(q);
Glocalnum = cellfun(@full, Glocalnum, 'UniformOutput',false);

parent = ax;
for i = 1:length(q)
    
    transforms{i} = hgtransform('parent', parent);
    parentPrev = parent;
    parent = transforms{i};
    
    transforms{i}.Matrix = Glocalnum{i};
    Glink1 = G0{i};
    P1 = Glink1(1:3, 4);
    Glink2 = Glink1*G0{i+1};
    P2 = Glink2(1:3, 4);
    
    G = G*Glocalnum{i};
    
    % geometric primitives
    if strcmpi(Jtype(i) , 'R')
        [jointX, jointY, jointZ] = createCylinder(options.jointradius, options.jointlength/2, options.jointlength/2, X(4:6, i), 0);
        % plot joint
        surf(jointX, jointY, jointZ, 'facecolor', 'r', 'edgecolor', 'none', 'parent', transforms{i});
        fill3(jointX(1,:), jointY(1,:),  jointZ(1,:), 'r', 'parent', transforms{i});
        fill3(jointX(2,:), jointY(2,:),  jointZ(2,:), 'r', 'parent', transforms{i});
        
        [linkX, linkY, linkZ] = createCylinder(options.jointradius/2, 0, -norm(P2-P1), (P2 - P1)./norm(P2-P1), 0);
    elseif strcmpi(Jtype(i) , 'P')
        [ptsP, faces] = createParallelepiped(options.jointradius, [options.jointlength*1.5 options.jointlength*1.5], X(1:3, i), 0);
        patch('faces', faces, 'vertices', ptsP.',  'facecolor', 'green', 'Parent', parentPrev);
        
        [linkX, linkY, linkZ] = createCylinder(options.jointradius/2, norm(P2-P1)-0.5, -norm(P2-P1), (P2 - P1)./norm(P2-P1), 0);
    end
    
    %plot link
    surf(linkX, linkY, linkZ, 'facecolor', 'b', 'edgecolor', 'none', 'parent', transforms{i});
    fill3(linkX(1,:), linkY(1,:),  linkZ(1,:), 'b', 'parent', transforms{i});
    fill3(linkX(2,:), linkY(2,:),  linkZ(2,:), 'b', 'parent', transforms{i});
end

% plot end effector
[EEX, EEY, EEZ] = createCylinder(options.jointradius/3, 0.08, 0, [0;0;1], 0);
G = G0{end};
EEX = G(1,1).*EEX + G(1,2).*EEY + G(1,3).*EEZ + G(1,4);
EEY = G(2,1).*EEX + G(2,2).*EEY + G(2,3).*EEZ + G(2,4);
EEZ = G(3,1).*EEX + G(3,2).*EEY + G(3,3).*EEZ + G(3,4);
surf(EEX, EEY, EEZ, 'facecolor', '#4DBEEE', 'edgecolor', 'none', 'parent', transforms{i});
patch(EEX(1,:), EEY(1,:), EEZ(1,:),'facecolor', '#4DBEEE', 'parent', transforms{i});
patch(EEX(2,:), EEY(2,:), EEZ(2,:),'facecolor', '#4DBEEE', 'parent', transforms{i});

end