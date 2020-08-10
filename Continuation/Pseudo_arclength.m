% This code does rigorous pseudo-arclength continuation in two directions,
% starting at a known periodic orbit (need to be enterred). The user is
% prompted to enter if it wants to save the curve as a csv file. If so, the
% path is enterred and the name of the .csv file corresponds to the first
% two coordinates (t0, x0) of the periodic orbit (same for figure). The
% periodic orbits in the .csv file are stored as column vectors and use the 
% normal ordering. The epsilon value is saved (4k+1 column) and the
% existence interval of the smooth curve is saved (next two columns). Since 
% continuation is done in two directions, the existence interval at a given 
% column corresponds to the curve between the periodic orbit at this column
% and the periodic orbit at the previous column, where previous means in
% the direction of the starting point. Here is an example:

% .csv file:
% periodic_orbit1 | epsilon_value1 | rmin1 | rmax1 |
% periodic orbit2 | epsilon_value2 | rmin2 | rmax2 |
% Initial_solution| Initial_epsilon|   0   |   0   |
% periodic_orbit3 | epsilon_value3 | rmin3 | rmax3 |
% periodic_orbit4 | epsilon_value4 | rmin4 | rmax4 |

% where rmin1 is between the 1st and second periodic orbit, rmin2 is
% between the initial solution(periodic orbit) and the 2nd, rmin3 is
% between the initial and 3rd solution and rmin4 is between the 3rd and 4th
% solution.

% Setting a value for "a" i.e. a positive, nonzero value (steepness of surface).
a = input('a = ');
% Setting a valid value for "e" i.e. a value in the interval [0, 1]
epsilon = input('e = ');
% Setting a valid value for "E" i.e. a nonnegative value (amplitude of oscillation).
E = input('E = ');
% Setting preferences for the search of periodic orbits
k = input('Enter an integer value for the number of impacts (bigger than 1): \n k = ');
p = input('Enter an integer value for the period: \n p = ');
save_curve_as_figure = input('Do you want to save the solution curve as a MATLAB .fig file? ');
if save_curve_as_figure == 1
    path_to_figs = input("Enter a path in single quotes (e.g. 'C:\Users\clmnt\documents') ");
end
save_curve_as_csv = input('Do you want to save the solution curve in a .csv file? ');
if save_curve_as_csv == 1
    path_to_csv = input('Enter a path in single quotes: ');
end
% Creating empty double array for plotting
em = double.empty;
% Edit the either of the following two lines (first one if solution is
% found in a .csv file. It must be saved as a column vector respecting the
% order (t0, x0, v0, w0, t1, x1, v1, w1, ..., tk-1, xk-1, vk-1, wk-1). The
% second line is if you want to explicitly enter the periodic orbit, again
% as a column vector respecting the order aforementioned.

solution_set = readmatrix('C:\Users\clmnt\Documents\McGill\Summer_2020\MATLAB\Faster_Billiard_2D\solutions\(1,1,1,2,2).csv');
%solution_set = [periodic orbit (column vector)]

% Do you want to refine the solution(s)?
Refining = input('Do you want to refine the inital solution? (Y=1/N=0) ');
% Number of times the parameter delta_s changes for 1 direction 
n = 3;
for s=1:size(solution_set,1)
    % Initial solution as a column vector
    X0 = [solution_set(s,1:4*k),epsilon].'
    % Refining solution if needed
    if Refining == 1
        f = @(X) zero_finding_problem(X(1:4*k),a,e,X(4*k+1),k,p);
        Df = @(X) jacobian(X(1:4*k),a,e,X(4*k+1),k,p);
        for i0=1:10
            X0(1:4*k) = X0(1:4*k) - Df(X0)\eye(4*k)*f(X0);
            if norm(f(X0),inf)< 0.5e-13
                break
            end
        end
    end
    if norm(f(X0),inf) >= 0.5e-13 || isnan(norm(f(X0),inf))
        error("The initial vector is not a solution.")
    end
    X_initial = X0;
    Xp = cell(n,2);
    Xm = cell(n,2);
    % Doing Rigorous parameter continuation in positive direction
    delta_s = 1e-1;
    for ip=1:n
        if ip>1
            if size(cell2mat(Xp(ip-1,1)),1) >= 2500
                delta_s = delta_s*2;
            else
                delta_s = delta_s/2;
            end
        end
        disp(['Iteration #', num2str(ip), ', the delta s parameter is now ', num2str(delta_s)])
        [solution,I,X0,pitchforkp,flag_backp] = Continuation(X0,delta_s,a,e,k,p);
        Xp(ip,1) = {solution};
        Xp(ip,2) = {I};
        if pitchforkp == 1 || flag_backp == 1
            break
        end       
    end 
    X0 = X_initial;
    % Doing Rigorous parameter continuation in negative direction
    delta_s = -1e-1;
    for im=1:n 
        if im>1
            if size(cell2mat(Xm(im-1,1)),1) >= 2500
                delta_s = delta_s*2;
            else
                delta_s = delta_s/2;
            end
        end
        disp(['Iteration #', num2str(im), ', the delta s parameter is now ', num2str(delta_s)])
        [solution,I,X0,pitchforkm,flag_backm] = Continuation(X0,delta_s,a,e,k,p);
        Xm(im,1) = {solution};
        Xm(im,2) = {I};
        if pitchforkm == 1 || flag_backm == 1
            break
        end
    end
    % Storing in an array
    X = em;
    for i=1:n
        X = [flip(cell2mat(Xm(i,1)),1),flip(cell2mat(Xm(i,2)),1); X];
        X = [X; cell2mat(Xp(i,1)),cell2mat(Xp(i,2))];
    end
    % Insert initial condition
    XM = cell2mat(Xm(:,1));
    X = [X(1:size(XM,1),:); [X_initial.',0,0]; X(size(XM,1)+1:end,:)];
    % Plotting result of both directions and saving figure
    name = [num2str(X_initial(1)),' and ',num2str(X_initial(2))];
    H = Plot_branches('E',a,e,epsilon,p,X(:,1:4*k+1));
    if save_curve_as_figure == 1
       savefig(H, [path_to_figs,name,'.fig']); 
    end
    % Saving big array with intervals
    if save_curve_as_csv == 1
        writematrix(X, [path_to_csv,name,'.csv']);
    end
    pause(3) 
    close(H)
end