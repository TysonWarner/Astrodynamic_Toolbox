function [v1,v2] = Lambert1(r1,r2,p,mu)
    r1n = norm(r1); r2n = norm(r2);
    %% Finding difference in true anomaly
    df = acos(dot(r1,r2)/(r1n*r2n));
    %% F&G Solutions
    % Finding velocity at position 1
    F = 1-(r2n/p)*(1-cos(df));
    G = ((r1n*r2n)/(sqrt(mu*p)))*sin(df);
    v1 = (r2-F.*r1)./G;
    % Finding velocity at position 2
    Fdot = ((dot(r1,v1))/(p*r1n))*(1-cos(df)) - sqrt(mu/(p*r1n^2))*sin(df);
    Gdot = 1 - (r1n/p)*(1-cos(df));
    v2 = Fdot.*r1 + Gdot.*v1;
end