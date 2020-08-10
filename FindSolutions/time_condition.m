function [y] = time_condition(periodic_orbit, a, e, E, k, p)
y(1) = true; % assume there are no zero in the interval
% The [y] array will contain the true/false values of each time segment
% in the solution (if respects time condition, true, else false)
% close the loop i.e. (tk, xk, vk, wk) = (t0+p, x0, v0, w0)
periodic_orbit(1+4*k:4+4*k) = [periodic_orbit(1)+p; periodic_orbit(2); ...
                                periodic_orbit(3); periodic_orbit(4)];    
for i=1:k % loop through the time components of the solution 
    if round(periodic_orbit(1+4*(i-1)),5)>=round(periodic_orbit(1+4*i),5) % if tn >= tn+1, solution is not physical
        y(i) = false;
        break
    end 
    % We make Tn a function of t (tn+1) only and input the rest as double
    f = @(t) next_time(t, periodic_orbit(1+4*(i-1)), periodic_orbit(2+4*(i-1)), ...
                        periodic_orbit(3+4*(i-1)), periodic_orbit(4+4*(i-1)), a, e, E);
    % (tn, tn+1)                
    interval = infsup(periodic_orbit(1+4*(i-1)), periodic_orbit(1+4*i));   
    % function from INTlab
    roots_function = verifynlssall(f, interval); % gives interval    
    if isempty(roots_function) ~= 1 % if roots array is not empty
        for j=1:length(roots_function) 
            x = true(length(roots_function), 1); % create a logical vector
            if round(periodic_orbit(1+4*(i-1)),3) <= round(sup(roots_function(j)),3) && round(periodic_orbit(1+4*(i-1)),3) >= round(inf(roots_function(j)),3)
                x(j) = 0; % take out the roots at infimum and supremum
            elseif round(periodic_orbit(1+4*i),3) <= round(sup(roots_function(j)),3) && round(periodic_orbit(1+4*i),3) >= round(inf(roots_function(j)),3)
                x(j) = false; % take out the roots at infimum and supremum            
            end
        end
    roots_function = roots_function(x); % take roots which are not at the boundaries    
    end
    if isempty(roots_function) ~= 1
        y(i) = false;
    end 
end
y = all(y);
end