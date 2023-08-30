function Jb = BodyJac(gst0, X, q)
%
% Jb = BodyJac(gst0,{Y1,var1},{Y2,var2},...,{Yn,varn});
% Funzione che calcola lo Jacobiano Body di un seriale
% a partire dalla parametrizzazione GLOBAL POE
%
%
% varargin : {Yn,varn} inserire tante celle 1x2 quanti sono i
% giunti del seriale;
% Yn : n esimo twist unitario del n esimo giunto (sdr spatial);
% varn : n esima variable di giunto;
%
% gst0: configurazione di riferimento iniziale
%
n = length(q);
cl = class(q);
g = gst0;
Jb = zeros(6, n, cl);

for i = n : -1 : 1
    
    g = expTw(X(:, i), q(i))*g;  %g_i+1,n
    Jb(:, i) = adjoint(rigidInverse(g))*X(:, i); % adjoint(g_i,n^-1)*Y_i
    
end

end