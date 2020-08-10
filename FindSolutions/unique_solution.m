function y = unique_solution(solution_set, periodic_orbit, k, p, PS)
% This function checks if the solution found has already been found
kp_array = cell2mat(solution_set(1:PS, 2)); % (k,p) of previous solutions
kp = [k, p]; % (k,p) of new solution
x = zeros(PS, 1); % assume (k,p) are all different
for i=1:PS
    x(i) = isequal(kp, kp_array(i,:)); % check for same (k,p) as the new solution 
end
x = logical(x); % logical array
y = ones(sum(x), 1); % assume the solution is different from all the other ones
if any(x)
    same_kp_solutions = round(cell2mat(solution_set(x,1).'), 4);
    for i=1:sum(x) % number of solutions for the same (k,p)
    y(i) = ~isequal(round(periodic_orbit, 4), same_kp_solutions(:,i)); % check if same solution
    end
end
y = all(y); % if unique then true, else false
end
