function [a,ecc,f,i,RAAN,AOP] = rv2oe(rE, vE, mu)
   
    %---------- Orbit Specific Energy ----------

    energy = ((norm(vE)^2)/2)-mu/(norm(rE));

    %---------- Semi-Major Axis ----------
    
    a = -0.5*mu/energy;
    
    %---------- Angular Momentum ----------
    
    h = cross(rE,vE);

    %---------- Eccentricity ----------
    
    e = (cross(vE,h)/mu)-(rE/norm(rE));
    ecc = norm(e);

    %---------- Inclination ----------
    
    i = acos((h(3))/norm(h));

    %---------- True Anomaly ----------
    
    if dot(rE,vE)>0
        f = acos((a*(1-ecc^2))/(norm(rE)*ecc) - 1/ecc);
    else
        f = -acos((a*(1-ecc^2))/(norm(rE)*ecc) - 1/ecc);
    end

    %---------- Right Ascension of Ascending Node ----------
    
    n = cross([0 0 1],h)/norm(cross([0 0 1],h));    % Node Vector
    sinRAAN = dot(n,[0 1 0]);
    cosRAAN = dot(n,[1 0 0]);
    if cosRAAN<0                                   % Finding RAAN
        RAAN = atan(sinRAAN/cosRAAN)+pi;
    else
        RAAN = atan(sinRAAN/cosRAAN);
    end
    
    %---------- Argument of Perigee ----------
    
    if dot(e,[0 0 1])<0
        AOP = -acos(dot(e,n)/norm(e));
    else
        AOP = acos(dot(e,n)/norm(e));
    end
end