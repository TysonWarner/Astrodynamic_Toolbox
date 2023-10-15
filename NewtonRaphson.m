function EccAnom = NewtonRaphson(a,e,mu,dt)

    delta = 0.0001;    % In rad
    
    % Calling function
    % E_new = NRI(a,e,u,dt,delta);
    
    % Getting final E_new
    % E_new_final = E_new(:,size(E_new,2));
    
    % ------------ Newton Raphson Iteration (NRI) -----------
    
    % Calculating M
    M = sqrt(mu/(a^3))*(dt);

    %Initializing
    error = 2 * delta;
    E(1) = M;
    i=1;

    %Looping to find E (or E_new)
    while error >= delta
        E(i+1) = E(i)-(M+e*sin(E(i))-E(i))/(e*cos(E(i))-1);
        error = abs(E(i+1) - E(i));
        i=i+1;
    end

    EccAnom = E(end);
end