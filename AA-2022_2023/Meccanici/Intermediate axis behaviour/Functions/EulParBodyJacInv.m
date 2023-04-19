function Jb_inv = EulParBodyJacInv(q)
%
% Jb_inv = EulParBodyJacInv(q)
%
% returns the inverse of the body jacobian corresponding to the Euler
% parameters (quaternion coefficients)
%
% input q: 4x1 array with the Euler parameters = [q0;q1;q2;q3]
%

q0 = q(1);
qvec = q(2:4);

Jb_inv = 0.5.*[-qvec.';...
    q0.*eye(3) + hat(qvec)];

end
