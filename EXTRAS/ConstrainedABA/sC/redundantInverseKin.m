function [Jr, b] = redundantInverseKin(J, w, q0_dot)
%
% [Jr, b] = redundantInverseKin(J, W, q0_dot)
%
% function that computes the weighted pseudo inverse jacobian (Jr) and the
% projection vector (b) of the joint velocities in the null of J 
%
% INPUTS: - J redundant robot jacobian
%         - w weight vector (weights for q_dot, we chose which joint velocities we want to minimize more)
%         - gradH objective function gradient to follow a reference (desired)  trajectory
%           defined by q0_dot
%         - k0 scaling coefficient for the objective H

Winv = diag(1./w);
JT = J.';
Jr = Winv*JT*inv(J*Winv*JT);
JrJ = Jr*J;
P = eye(size(JrJ, 1)) - JrJ;
b = P*q0_dot; % null projection of q0_dot

end