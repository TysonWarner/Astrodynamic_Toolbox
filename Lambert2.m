function [dt] = Lambert2(r1,v1,df,mu)
r1n = norm(r1); v1n = norm(r2);
%% Orbital Elements
energy = (v1n^2)/2 - mu/r;
a = -mu/(2*energy);
h = norm(cross(r1,v1));
p = (h^2)/mu;
if energy<0
    orbit = 'elliptic';
else
    orbit = 'hyperbolic';
end
a = -mu/(2*energy);
h = norm(cross(r1,v1));
p = (h^2)/mu;
%% F&G Solutions

if strcmp(orbit,'elliptic') == 1

else

end
% Finding velocity at position 1
F = 1-(r2n/p)*(1-cos(df));
G = ((r1n*r2n)/(sqrt(mu*p)))*sin(df);
v1 = (r2-F.*r1)./G;
% Finding velocity at position 2
Fdot = ((dot(r1,v1))/(p*r1n))*(1-cos(df)) - sqrt(mu/(p*r1n^2))*sin(df);
Gdot = 1 - (r1n/p)*(1-cos(df));
v2 = Fdot.*r1 + Gdot.*v1;
end