% Setting a value for "a" i.e. a positive, nonzero value (steepness of surface).
a = input('a = ');
% Setting a valid value for "e" i.e. a value in the interval [0, 1]
e = input('e = ');
% Setting a valid value for "E" i.e. a nonnegative value (amplitude of oscillation).
E = input('E = ');
% Setting preferences for the search of periodic orbits
set_k_manually = input('Do you want to fix k (Yes = 1/No = 0)?');
if set_k_manually == 1
    k = input('Enter an integer value (bigger than 1): \n k = ');
end
set_p_manually = input('Do you want to fix p (Yes = 1/No = 0)?');
if set_p_manually == 1
    p = input('Enter an integer value: \n p = ');
end
n_PS = input('Number of solutions to find: ');
% This will be useful since the periodic solutions can be of different size
% 1st column: coordinates; 2nd column: (k,p); 3rd column: EI
solution_set = cell(n_PS, 3);
% Number of periodic solutions found
PS = 0;
disp(['Searching for ', num2str(n_PS), ' periodic orbit(s)...'])
while PS< n_PS
    if set_p_manually ~= 1
        p = randi(10); % 1<=p<=10
    end
    if set_k_manually ~= 1 
        k = randi(9)+1; % 2<=k<=10
    end
    periodic_orbit = ones(4*k, 1); % preallocating array size
    % Generating (possible) random initial conditions
    periodic_orbit(1) = rand(1);
    periodic_orbit(2) = -5+10*rand(1); % -5 <= x0 <= 5
    periodic_orbit(3) = sign(periodic_orbit(2))*(5*rand(1)); % |v0| <= 5
    periodic_orbit(4) = -20+40*rand(1); % -20 <= w0 <= 20
    for i=1:k-1
        periodic_orbit(1+4*i) = periodic_orbit(1+4*(i-1)) + i*p/k;
        periodic_orbit(2+4*i) = -5+10*rand(1);
        periodic_orbit(3+4*i) = sign(periodic_orbit(2+4*i))*(5*rand(1));
        periodic_orbit(4+4*i) = -10+20*rand(1);
    end
    for i=1:100 % iterating with Newton's operator
        % Newton's method
        periodic_orbit = periodic_orbit - jacobian(periodic_orbit,a,e,E,k,p)\eye(4*k)*zero_finding_problem(periodic_orbit,a,e,E,k,p);
        % if smaller than given tolerance (tol = 1e-12) 
        if norm(zero_finding_problem(periodic_orbit,a,e,E,k,p),inf)<=1e-13
            % readjust time components such that 0 <= t0 <= 1
            periodic_orbit = normalized_t0(periodic_orbit, k);
            % if it satsifies the time ordering AND minimum time
            if time_condition(periodic_orbit,a,e,E,k,p) == 0
                break
            else % if it satisfies minimum time
                if unique_solution(solution_set, periodic_orbit, k, p, PS) % if solution is unique
                    PS = PS + 1;
                    disp(['n = ', num2str(PS), '. A periodic orbit has been found.'])
                    solution_set(PS, 1) = {periodic_orbit};
                    solution_set(PS, 2) = {[k,p]};
                    % Store the existence interval in third column
                    disp('Proving the solution.')
                    solution_set(PS, 3) = {existence_interval(periodic_orbit,a,e,E,k,p)};
                    fprintf(['Searching for ', num2str(n_PS-PS), ' other periodic orbit(s)...\n\n'])
                    break
                else % solution is not unique, generate new initial conditions
                    break
                end
            end
        % if F is diverging, stop and generate new initial conditions
        elseif norm(zero_finding_problem(periodic_orbit,a,e,E,k,p),inf) >1e5 || i==100
            break
        end          
    end
end  
disp('Search completed.')
% Creating a table with the solutions, the (k,p) values and their EI
T = cell2table(solution_set);   
save2csv = input('Save to csv file? (Yes = 1/No = 0)');
if save2csv == 1
    name = input("Write the name of the csv file e.g. 'name.csv'):");
    writetable(T, name);
end