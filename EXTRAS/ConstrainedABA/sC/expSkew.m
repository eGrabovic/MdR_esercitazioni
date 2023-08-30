function X = expSkew(axis, th)
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

absAxis = abs(axis);

if absAxis(1) <= 1e-8 && absAxis(2) <= 1e-8
    
    if axis(3) > 0
        
        X = rotZ(th);
        
    elseif axis(3) < 0
        
        X = rotZ(-th);
    end
    
elseif absAxis(2) <= 1e-8 && absAxis(3) <= 1e-8
    
    if axis(1) > 0
        
        X = rotX(th);
        
    elseif axis(1) < 0
        
        X = rotX(-th);
    end
    
elseif absAxis(1) <= 1e-8 && absAxis(3) <= 1e-8
    
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
