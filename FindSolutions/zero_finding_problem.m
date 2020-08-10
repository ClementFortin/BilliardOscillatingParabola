function [F] = zero_finding_problem(periodic_orbit, a, e, E, k, p)
% This function corresponds to the "f" map
% periodic_orbit = (t0, x0, v0, w0, ..., tk, xk, vk, wk)
% close the loop i.e. (tk, xk, vk, wk) = (t0+p, x0, v0, w0)
periodic_orbit(1+4*k:4+4*k) = [periodic_orbit(1)+p; periodic_orbit(2); ...
                                periodic_orbit(3); periodic_orbit(4)];
for i=1:k
    % Outputs of the reset law
    R_out(1+2*(i-1):2+2*(i-1)) = R_map(periodic_orbit(1+4*(i-1)), ...
                    periodic_orbit(2+4*(i-1)), periodic_orbit(3+4*(i-1)), ...
                    periodic_orbit(4+4*(i-1)), a, E, e);
    % Outputs of the S map
    S_out(1+4*(i-1):4+4*(i-1)) = S_map(periodic_orbit(1+4*i), ...
                    periodic_orbit(1+4*(i-1)), periodic_orbit(2+4*(i-1)), ...
                    R_out(1+2*(i-1)), R_out(2+2*(i-1)), a, E);
    % Zero-finding problem
    F(1+4*(i-1):4+4*(i-1)) = [next_time(periodic_orbit(1+4*i), periodic_orbit(1+4*(i-1)), ...
                              periodic_orbit(2+4*(i-1)), periodic_orbit(3+4*(i-1)), ...
                              periodic_orbit(4+4*(i-1)), a, e, E); % the appended zero-finding problem
                              S_out(1+4*(i-1)) - periodic_orbit(2+4*i); % S(xn) - xn+1
                              S_out(3+4*(i-1)) - periodic_orbit(3+4*i); % S(vn) - vn+1
                              S_out(4+4*(i-1)) - periodic_orbit(4+4*i)]; %S(wn) - wn+1

end
F = F(:);
end