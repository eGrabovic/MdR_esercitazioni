function out = graspMatrix(p)

I = eye(3);
out = [     I,  zeros(3,3);...
       hat(p),          I];

end