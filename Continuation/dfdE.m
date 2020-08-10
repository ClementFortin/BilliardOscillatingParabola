function df = dfdE(X,a,e,k,p)
% This function computes the partial derivatives of the zero-finding
% problem w.r.t to the epsilon parameter (E)
E = X(1+4*k);
X(1+4*k:4+4*k) = [X(1)+p; X(2); X(3); X(4)];
X(1+4*(k+1)) = E;
% Partials w.r.t to the output of the R_map
dv = @(n) 2*pi*cos(2*pi*X(1+4*n))*(e+1)*(-2*a*X(2+4*n))/(1+4*a^2*X(2+4*n)^2);
dw = @(n) 2*pi*cos(2*pi*X(1+4*n))*(e+1)/(1+4*a^2*X(2+4*n)^2);

% Ouputs of the R_map
for n=0:k-1
    
    R_out = R_map(X(1+4*n),X(2+4*n),X(3+4*n),X(4+4*n),a,X(1+4*(k+1)),e);

    df(1+4*n:4*(n+1)) = [2*a*(X(2+4*n)+R_out(1)*(X(1+4*(n+1))-X(1+4*n)))*dv(n)*(X(1+4*(n+1))-X(1+4*n))+sin(2*pi*X(1+4*(n+1)))-sin(2*pi*X(1+4*n))-dw(n)*(X(1+4*(n+1))-X(1+4*n));
                        dv(n)*(X(1+4*(n+1))-X(1+4*n));
                        dv(n);
                        dw(n)
                        ];
end
df = df(:);
end