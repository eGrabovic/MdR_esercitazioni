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