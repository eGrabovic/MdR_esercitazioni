function dJ_dq = CGDHjacobianDerivative(T0j, Tj_dq, CG_i, bodyIndex)
%
%
%


cl = class(T0j{1});
Jac = zeros(6, length(T0j), cl);
Jac_der = cell(1, length(T0j), cl);
n = body_index;


end