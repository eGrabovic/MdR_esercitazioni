function out = rodriguezAxisAngle(ax, x)
%
% Rodriguez formula for the rotation matrix with axis angle
% parametrization
%

h_ax = hat(ax);
out = (eye(3) + h_ax.*sin(x) + h_ax*h_ax.*(1 - cos(x)));

end