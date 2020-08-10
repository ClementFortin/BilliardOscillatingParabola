function [X0,I,Last,pitchfork,flag_back] = Continuation(X0,initial_delta_s,a,e,k,p)
% This function does continuation in the 'E' parameter
X_initial = X0;
pitchfork = 0;
% -------------------------------------------------
% ZERO-FINDING PROBLEM AND jacobian MATRIX
% -------------------------------------------------
% X is a 4k+1 x 1 vector
f = @(X) zero_finding_problem(X(1:4*k),a,e,X(4*k+1),k,p);
Df = @(X) jacobian(X(1:4*k),a,e,X(4*k+1),k,p);
% --------------------------------------------------
% UNIT TANGENT VECTOR OBTAINED VIA NULL FUNCTION
% --------------------------------------------------
% We want to find a vector in the null space of the jacobian
Df_4kx4k1 = @(X) [Df(X),dfdE(X,a,e,k,p)];
% 4k+1 dimension tangent vector 
X_tangent = @(X) null(Df_4kx4k1(X));
% ---------------------------------------------------
% CONSTRUCTING THE NEWTON OPERATOR
% ---------------------------------------------------
% Predictor
X_predictor = @(X0,X0_dot,delta_s) X0 + delta_s.*X0_dot;
plane = @(X,X0,X0_dot,delta_s) (X-X_predictor(X0,X0_dot,delta_s)).'*X0_dot;
% New zero-finding problem with the plane
F = @(X,X0,X0_dot,delta_s) [plane(X,X0,X0_dot,delta_s);f(X)];
% Constructing the jacobian matrix with the plane
DF = @(X,X0_dot) [X0_dot.'; Df_4kx4k1(X)];
% ---------------------------------------------------
% PARAMETER CONTINUATION (PSEUDO-ARCLENGTH)
% ---------------------------------------------------
% number of steps for parameter continuation
n = 5000;
% Changing array to type double
X = zeros(4*k+1,n);
X_dot = zeros(4*k+1,n); 
delta_s = zeros(1,n);
% Obtaining initial tangent vector 
X_dot(:,1) = X_tangent(X0(:,1));
% Initial step size
delta_s(1) = initial_delta_s;
% Declaring array for existence interval of the "corrected" solution
I = zeros(n,2);
rmin = zeros(1,n);
rmax = zeros(1,n);
rho = zeros(2,n);
% flag to break out of for loops
flag = 0; 
flag_back = 0;
% Newton's method
for j=1:n
    % Initializing with predictor
    X(:,1) = X_predictor(X0(:,j),X_dot(:,j),delta_s(j));
    % Newton's method
    for i=1:100
        X(:,i+1) = X(:,i) - (DF(X(:,i),X_dot(:,j))\eye(4*k+1))*F(X(:,i),X0(:,j),X_dot(:,j),delta_s(j));
        % If tolerance is attained start proving
        if norm(F(X(:,i),X0(:,j),X_dot(:,j),delta_s(j)),inf)<=0.5e-13
             % ---------------------------------------------------
             % SETTING NEW VALUES
             % ---------------------------------------------------
             % Storing the value as a solution
             X0(:,j+1) = X(:,i);
             if j>1 % If the new solution is in between the two previous solution for all indices, stop
                 for idx=1:size(X0(:,j+1))
                     if round(X0(idx,j+1),13) <= round(max(X0(idx,j),X0(idx,j-1)),13) && round(X0(idx,j+1),13) >= round(min(X0(idx,j),X0(idx,j-1)),13)
                         m(idx) = 1;
                     else 
                         m(idx) = 0;
                     end
                 end
                 if all(m)
                     flag_back = 1;
                     break
                 end
             end
             % Storing the tangent vector at that solution (same 
             % direction as previous one). If we are at a pitchfork, 
             % take the tangent vector which is "more" parallel to the
             % previous one
             if size(X_tangent(X0(:,j+1)),2) ~= 1
                disp(['There is a pitchfork bifurcation nearby.',newline])
                pitchfork = 1;
                tangent_vectors = X_tangent(X0(:,j+1));
                for i3=1:size(X_tangent(X0(:,j+1)),2)
                   dot_product(i3) = dot(tangent_vectors(:,i3),X_dot(:,j));
                end
                [~,index] = max(abs(dot_product));
                X_dot(:,j+1) = sign(X_dot(:,j).'*tangent_vectors(:,index))*tangent_vectors(:,index);
            else
                X_dot(:,j+1) = sign(X_dot(:,j).'*X_tangent(X0(:,j+1)))*X_tangent(X0(:,j+1));
            end
            % Storing the number of steps for new delta_s
            delta_s(j+1) = delta_s(j)*2^((4-i)/3);
            % Setting a maximum step size 
            if abs(delta_s(j+1)) > abs(delta_s(1))
                delta_s(j+1) = delta_s(1);
            end
            % ---------------------------------------------------
            % PROVING EXISTENCE OF SMOOTH CURVE
            % ---------------------------------------------------
            [I(j+1,:),flag,~] = solution_curve(X0(:,j),X0(:,j+1),X_dot(:,j),X_dot(:,j+1),a,e,k,p,j);
            break
        elseif norm(F(X(:,i),X0(:,j),X_dot(:,j),delta_s(j)),inf)>1e3 || i==100 || flag == 1
            disp(['The sequence did NOT converge for j = ', num2str(j)])
            flag = 1;
            break
        end               
    end
    if flag == 1 || flag_back == 1
        break
    end
end
if j~=n % Truncate arrays
   X0(:,j+1:end) = [];
   X0(:,1) = [];
   I(j+1:end,:) = [];
   I(1,:) = [];
end
if isempty(X0)
    Last = X_initial;
else
    Last = X0(:,end);
end
X0 = X0.';
end


function [I, flag, X1_dot] = solution_curve(X0, X1, X0_dot, X1_dot,a, e, k, p, i)
    % This function attempts to prove that there exists a smooth branch curve
    % between the solution X0 and the solution X1 obtained via pseudo-arclength
    % continuation.
    % Setting the balls radius in which we want to prove existence
    r_star = 1e-4;
    s = infsup(0,1);
    f = @(X) zero_finding_problem(X(1:4*k),a,e,X(4*k+1),k,p);
    Df = @(X) jacobian(X(1:4*k),a,e,X(4*k+1),k,p);
    % We want to find a vector in the null space of this matrix
    Df_4kx4k1 = @(X) [Df(X),dfdE(X,a,e,k,p)];
    % 4k+1 dimension tangent vector, if not specified already
    if isempty(X1_dot)
        X_tangent = @(X) null(Df_4kx4k1(X));
        % If at a pitchfork bifurcation, take the MOST parallel tangent vector
        if size(X_tangent(X1),2) ~= 1
            disp('There is a pitchfork bifurcation nearby.')   
            tangent_vectors = X_tangent(X1);
            for j=1:size(X_tangent(X1),2)
                dot_product(j) = dot(tangent_vectors(:,j),X0_dot);
            end
            [~,index] = max(abs(dot_product));
            X1_dot = sign(X0_dot.'*tangent_vectors(:,index))*tangent_vectors(:,index);
        else
            X1_dot = sign(X0_dot.'*X_tangent(X1))*X_tangent(X1);
        end
    end
    % Constructing the jacobian matrix with the plane
    DF = @(X,X0_dot) [X0_dot.'; Df_4kx4k1(X)];
    DDF = @(X) D2F(X,a,e,p);

    Xs = X0 + s*(X1-X0);
    Xs_dot = X0_dot + s*(X1_dot-X0_dot);
    plane_0 = @(X) (X-X0).'*X0_dot; % s = 0
    F0 = @(X) [plane_0(X); f(X)]; % s = 0 

    % Linear operator
    A = DF(X0,X0_dot)\eye(4*k+1);
    % Y0 bound
    Y0 = norm(A*F0(X0), inf);
    % Use fundamental theorem of calculus (integral) to set Y0_hat bound
    integrand = @(s) [zeros(1,4*k+1); Df_4kx4k1(X0+s*(X1-X0))];
    DF_integrated = integral(integrand, 0, 1, 'ArrayValued', true);
    Y0_hat = sup(norm(A*(DF_integrated*(X1-X0)), inf));
    % Set Z0 bound directly
    Z0 = norm(eye(4*k+1) - A*DF(X0,X0_dot),inf);
    % Use mean value theorem to set Z0_hat bound i.e. DFs(Xs) - DF(X0) = s[D2F(C)(X1-X0)] = s[DF(X1)-DF(X0)]
    First_row_diff_jacobians = (X1_dot-X0_dot).';
    Other_rows_diff_jacobians = Df_4kx4k1(X1)-Df_4kx4k1(X0);
    diff_jacobians = [First_row_diff_jacobians; Other_rows_diff_jacobians];
    Z0_hat = norm(A*(diff_jacobians),inf);
    % We set Z2 to be a local bound (using r_star and eq. 2.35)
    radius = rad(Xs) + r_star; 
    b = midrad(mid(Xs),radius);
    Z2 = sup(norm(A*(DDF(b)),inf));
    % roots of the polynomial
    polynomial = [Z2 -(1-Z0_hat-Z0) Y0+Y0_hat];
    rho = roots(polynomial);
    % verify that the roots are real positive numbers
    if all(imag(rho)==0) && min(rho) < r_star && max(rho) > 0 
        rmin = min(rho);
        rmax = min([max(rho), r_star]);
        I = [rmin, rmax];
        flag = 0;
    else
        % Use normal Z2 bound (w/ jacobian difference)
        Z2 = sup(norm(A*(DF(b,Xs_dot)-DF(Xs,Xs_dot)),inf));
        polynomial = [Z2 -(1-Z0_hat-Z0) Y0+Y0_hat];
        rho = roots(polynomial);
        % verify that the roots are real positive numbers
        if all(imag(rho)==0) && min(rho) < r_star && max(rho) > 0 
            rmin = min(rho);
            rmax = min([max(rho), r_star]);
            I = [rmin, rmax];
            flag = 0;
        else
            disp('The existence interval of the solution branch could not be defined.')
            flag = 1;
            I = [0,0];
        end
    end
    if flag == 0
        smoothness = smooth_curve(X0,X1,X0_dot,X1_dot,a,e,k,p,rmin);
        if smoothness == 1
            disp(['The pseudo-arclength exists and is smooth for i=', num2str(i), '.'])
        elseif smoothness ~= 1
            disp(['The pseudo-arclength exists but is not smooth for i=', num2str(i), '.'])
            flag = 1;
        end
    end
end


function smoothness = smooth_curve(X0, X1, X0_dot, X1_dot, a, e, k, p, rmin)
    % This function check the smoothness of the solution curve between X0 and X1
    if isempty(X1_dot)
        Df = @(X) jacobian(X(1:4*k),a,e,X(4*k+1),k,p);
        % We want to find a vector in the null space of this matrix
        Df_4kx4k1 = @(X) [Df(X),dfdE(X,a,e,k,p)];
        X_tangent = @(X) null(Df_4kx4k1(X));
        % If at a pitchfork bifurcation, take the MOST parallel tangent vector
        if size(X_tangent(X1),2) ~= 1
            disp('There is a pitchfork bifurcation nearby.')   
            tangent_vectors = X_tangent(X1);
            for j=1:size(X_tangent(X1),2)
                dot_product(j) = dot(tangent_vectors(:,j),X0_dot);
            end
            [~,index] = max(abs(dot_product));
            X1_dot = sign(X0_dot.'*tangent_vectors(:,index))*tangent_vectors(:,index);
        else
            X1_dot = sign(X0_dot.'*X_tangent(X1))*X_tangent(X1);
        end
    end
    % Defining the delta vector quantities
    delta_bar = X1 - X0;
    delta_dot = X1_dot-X0_dot;
    % Verifying that the path is regular (no secondary bifurcations)
    result = sign(-abs(delta_bar.'*X0_dot)+rmin*norm(delta_dot,inf)+abs(delta_bar.'*delta_dot));
    if result == 1 || result == 0
        smoothness = 0;
    elseif result == -1
        smoothness = 1;
    end
end
