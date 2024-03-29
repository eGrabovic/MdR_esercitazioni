function out = omegaFromQuat2(Q, Qdot)
R = quatToRot(Q);
q = Q(2:4);
q0 = Q(1);
Rdot = 2*hat(q)*Qdot(1) + (2*q0*eye(3) + 4*hat(q))*...
    (hat([1;0;0]).*Qdot(2) + hat([0;1;0]).*Qdot(3) + hat([0;0;1]).*Qdot(4));
out = vecForm(Rdot*R.');
end