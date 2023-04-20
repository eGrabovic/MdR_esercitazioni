function [Gglobal, Grel] = FWkin_localPOE(Goffset, X, q, options)
%
% inputs: Goffset: 1xn cell of 4x4 offset matrices 
%
%         X: 6xn matrix where each column is the unitary twist of each joint
%
%         q: nx1 array where each element is the joint value
%
% options: 'jstart' & 'jend': starting and ending joint for the forward
%                             kinematic computation
%
% outputs the forward kinematics with local POE formulation and each local
% relative transforms

arguments
    Goffset
    X
    q
    options.jstart = 1
    options.jend = length(q)
    options.EEoffset = eye(4);
end

jstart = options.jstart;
jend = options.jend;

Grel = cell(1, jend);
Gglobal = eye(4);

for i = jstart:jend
    Grel{i} = Goffset{i}*expTw(X(:, i), q(i));
    Gglobal = Gglobal*Grel{i};
end
Gglobal = Gglobal*options.EEoffset;
end