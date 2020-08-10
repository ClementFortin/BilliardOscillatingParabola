function y = normalized_t0(periodic_orbit, k)
% This function normalizes the time at which the motion starts (0<=t0<1)
% if t0 is not between 0 and 1
if periodic_orbit(1) < 0 || periodic_orbit(1) >= 1
    % then we must add this factor to every time values in periodic_orbit
    factor = -floor(periodic_orbit(1));
    for i=1:k
        periodic_orbit(1+4*(i-1)) = periodic_orbit(1+4*(i-1)) + factor;
    end
end
y = periodic_orbit; % return solution with adjusted time values
end