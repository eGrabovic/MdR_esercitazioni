function J = spatialJac(gst0, X, q)
%
% Jb = spatialJac(gst0,{Y1,var1},{Y2,var2},...,{Yn,varn});
% Funzione che calcola lo Jacobiano spatial di un seriale
%
%
% varargin : {Yn,varn} inserire tante celle 1x2 quanti sono i
% giunti del seriale;
% Yn : n esimo twist unitario del n esimo giunto;
% varn : n esima variable di configurazione del giunto;
%
% gst0: configurazione di riferimento iniziale
%
n = length(q);

g = expTw(X(:, 1), q(1));
J = zeros(6, n, class(q));
J(:, 1) = X(:, 1);
for i = 2:n
    J(:, i) = adjoint(g)*X(:, i);
    g = g*expTw(X(:, i), q(i));
end

end