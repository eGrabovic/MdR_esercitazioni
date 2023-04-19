function out = axisAngleToQuat(theta, n)
%
% out = axisAngleToQuat(theta, n)
%
% quaternion parameters from axis angle 
%
out = zeros(4,1, class(theta));
out(1) = cos(theta/2);
out(2:4) = sin(theta/2).*n;
end