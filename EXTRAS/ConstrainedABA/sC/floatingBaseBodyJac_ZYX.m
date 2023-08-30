function Jac = floatingBaseBodyJac_ZYX(q)

alpha = q(4:6);
cl = class(q);
Jw = bodyJac_ZYX(alpha);
R = rotZ(alpha(1))*rotY(alpha(2))*rotX(alpha(3));
O = zeros(3,3,cl);
Jac = [R.',   O;...
         O,  Jw];

end