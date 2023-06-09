function X = expSkew(axis,th)
%
% Calcola l'esponenziale e^(omega*theta);
% omega = versore asse (axis);
% theta = angolo (th);
%
% utilizza la formula di Rodriguez quando viene assegnato un
% asse generico;
%
% vengono direttamente riportate invece le rotazioni elementari
% nel caso l'asse coincida con uno degli assi principali (per
% velocizzare i calcoli)

if abs(axis.'*axis) -1 >= 1e-8
    error('axis requires to be a versor');
end

if abs(axis(1)) <= 1e-8 && abs(axis(2)) <= 1e-8
    
    if axis(3) > 0
        
        X = rotZ(th);
        
    elseif axis(3) < 0
        
        X = rotZ(-th);
    end
    
elseif abs(axis(2)) <= 1e-8 && abs(axis(3)) <= 1e-8
    
    if axis(1) > 0
        
        X = rotX(th);
        
    elseif axis(1) < 0
        
        X = rotX(-th);
    end
    
elseif abs(axis(1)) <= 1e-8 && abs(axis(3)) <= 1e-8
    
    if axis(2) > 0
        
        X = rotY(th);
        
    elseif axis(2) < 0
        
        X = rotY(-th);
    end
else
    axisHat = hat(axis);
    X = eye(3) + axisHat*sin(th) + axisHat*axisHat*(1 - cos(th));
    %rodriguez formula
end
