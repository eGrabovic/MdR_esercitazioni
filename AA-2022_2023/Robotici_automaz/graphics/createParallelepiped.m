function [points, faces] = createParallelepiped(edge, len, axis, axialOff)
%
% -edge = length of edges that define the cross-section rectangle (2
%   values), square otherwise
% -len = length of the parallelepiped (2 values) one for positive extension over
%      origin, the other for negative, otherwise intermediate to origin (if
%      1 value is given)
%
%      |    len(1)   |     len(2)  |
%      ____________________________
%      |             |             |
%  ------------------O-----------------> axis
%      |             |             |
%      ____________________________
%                    |
%                    |
%                    v 
a = edge;
if length(a) == 1
    a = repmat(a, 1, 2);
end
if length(len) == 2
    h1 = len(1);
    h2 = len(2);
else
    h1 = 0;
    h2 = len;
end
P1 = [a(1)/2;-a(2)/2;-h1]; P2 = [a(1)/2;a(2)/2;-h1]; P3 = [a(1)/2;a(2)/2;h2]; P4 = [a(1)/2;-a(2)/2;h2];
P5 = [-a(1)/2;a(2)/2;-h1]; P6 = [-a(1)/2;a(2)/2;h2]; P7 = [-a(1)/2;-a(2)/2;-h1]; P8 = [-a(1)/2;-a(2)/2;h2];
points = [P1,P2,P3,P4,P5,P6,P7,P8] - [0; 0; axialOff];
theta = acos(axis(3));
phi = atan2(axis(1), axis(2));
R = rotZ(phi)*rotX(theta);
points = R*points;
faces = [1,2,3,4;2,5,6,3; 3,6,8,4; 8,6,5,7; 1,4,8,7; 1,7,5,2];

end