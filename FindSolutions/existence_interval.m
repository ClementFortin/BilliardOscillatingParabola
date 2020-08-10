function I = existence_interval(periodic_orbit,a,e,E,k,p)
% A is the linear operator i.e. DF(numerical_solution)^-1
A = jacobian(periodic_orbit,a,e,E,k,p)\eye(4*k);
% we iterate from r = 1e-9 with step size 1e-5 to find Z just smaller than 1
% when Z>=1, we take the previous value of r such that Z<1 and iterate with
% step size 1e-7 and do the same with 1e-9 afterwards.
r = ones(100000,1);
r(1) = 1e-16;
for i=2:100000
    r(i) = 1e-9 + (i-1)*1e-5;
    ball_of_radius_r = midrad(periodic_orbit, r(i));
    DF = jacobian(ball_of_radius_r,a,e,E,k,p);
    DT_interval = eye(4*k) - A*DF;
    Z = sup(norm(DT_interval, inf));
    if Z>=1
        for j=i:100000
            r(j) = r(j-1) + j*1e-7;
            ball_of_radius_r = midrad(periodic_orbit, r(j));
            DF = jacobian(ball_of_radius_r,a,e,E,k,p);
            DT_interval = eye(4*k) - A*DF;
            Z = sup(norm(DT_interval, inf));
            if Z>=1
                for q=j:100000
                    r(q) = r(q-1) + q*1e-9;
                    ball_of_radius_r = midrad(periodic_orbit, r(q));
                    DF = jacobian(ball_of_radius_r,a,e,E,k,p);
                    DT_interval = eye(4*k) - A*DF;
                    Z = sup(norm(DT_interval, inf));
                    if Z>=1
                        break
                    end
                end
                break
            end
        end
        break
    end
end
r_star = r(q-1); % just before Z>=1
ball_of_radius_r_star = midrad(periodic_orbit, r_star);
DF = jacobian(ball_of_radius_r_star,a,e,E,k,p);
DT_interval = eye(4*k) - A*DF;
Z = sup(norm(DT_interval, inf));
% Z has been set. 
Y = norm(A*zero_finding_problem(periodic_orbit,a,e,E,k,p), inf);
rho = [Z-1 Y]; % the linear polynomial
root = roots(rho); % the root of the linear polynomial
I = [root, r_star]; % the existence interval
end