function T = next_time(t, s, x, v, w, a, e, E)
% This function is the "T" map
% periodic_orbit = (t0, x0, v0, w0, ..., tk, xk, vk, wk)
% Outputs of the reset law for the i-th impact
R_out = R_map(s, x, v, w, a, E, e);
% Zero-finding problem which tn+1 must satisfy
T = a.*((x + R_out(1).*(t-s)).^2 - x.^2) + E.*(sin(2.*pi.*t) - sin(2.*pi.*s)) ... 
        - R_out(2).*(t-s) + 4.905.*(t-s).^2;
end