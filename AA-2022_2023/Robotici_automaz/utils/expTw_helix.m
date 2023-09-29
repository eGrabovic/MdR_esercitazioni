function E = expTw_helix(unitTwist, th, h)
%
% Exponential twist computation
%
% unitTwist = unitary twist of the joint
%
% th = joint variable value
%


if h == Inf % prismatic
    
    axis = unitTwist(1:3);
    R = eye(3);
    d = axis*th;
    
elseif h == 0 % revolute
    
    axis = unitTwist(4:6);
    v = unitTwist(1:3);
    q = -cross(v,axis);
    R = expSkew(axis,th);
    d =(eye(3) - R)*q;
    
else
    
    axis = unitTwist(4:6);
    v = unitTwist(1:3);
    q = -cross(v,axis);
    R = expSkew(axis,th);
    d =(eye(3) - R)*q + axis.*h.*th;
    
end

E = [R d; 0 0 0 1];

end