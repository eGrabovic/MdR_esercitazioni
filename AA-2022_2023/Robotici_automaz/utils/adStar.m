function out = adStar(twist, varargin)
%
%
%
w = twist(4:6);
v = twist(1:3);
what = hat(w);
out = [what, zeros(3, 3);hat(v), what];

end