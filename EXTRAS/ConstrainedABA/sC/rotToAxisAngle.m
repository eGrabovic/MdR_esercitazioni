function [n, theta] = rotToAxisAngle(R)

theta = acos((trace(R)-1)/2) ;
n = 1/(2*sin(theta)+1e-8).*(vecForm(R)-vecForm(R')); 

end