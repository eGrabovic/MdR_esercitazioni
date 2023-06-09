function J = spatialJac_localPOE(Goffset, X, q, options)
%
% 
%
arguments
    Goffset
    X
    q
    options.EEoffset = eye(4);
    options.Toffset0 = eye(4);
end

n = length(q);
EEoffset = options.EEoffset;
Goffset{end} = Goffset{end}*EEoffset;
Goffset{1} = options.Toffset0*Goffset{1};
X(:, end) = adjointInv(EEoffset)*X(:, end);
J = zeros(6, n, class(q(1)));
B = eye(4);

for i = 1:1:n
    B = B*Goffset{i};
    J(:,i) = adjoint(B)*X(:, i);
    B = B*expTw(X(:, i), q(i));
end

end