function [xn_prime] = R_map(t, x, v, w, a, E, e)
% This function is the reset law
xn_prime = (1/(1+(2*a*x)^2))*[1-e*(2*a*x)^2, 2*a*x*(1+e); 2*a*x*(1+e), ...
    (2*a*x)^2-e]*[v; w] + (2*pi*E*cos(2*pi*t)*(1+e)/(1+(2*a*x)^2))*[-2*a*x; 1];
end



