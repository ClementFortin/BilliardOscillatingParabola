function y = jacobian(periodic_orbit, a, e, E, k, p)
% This function creates two matrices (block1, block2) for each n. Suppose
% k = 3, then Dim(Df) = 12x12. However, there are a lot of zero entries. 
% The first block contains the "n" partial derivatives associated with the
% "n-th" impact and the second contains the "n+1-th" partial derivatives
% associated with the "n-th" impact. The rest only contains zeros. For k=3,
% this looks like
%
% Df = 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%                       %                      %                       %
%         Block 1       %        Block 2       %            0          %
%                       %                      %                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       %                      %                       %
%            0          %        Block 1       %         Block 2       %
%                       %                      %                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       %                      %                       %
%         Block 2       %           0          %         Block 1       %
%                       %                      %                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% where each block is of dimension 4x4.
% (tk, xk, vk, wk) = (t0+p, x0, v0, w0)
periodic_orbit(1+4*k:4+4*k) = [periodic_orbit(1)+p; periodic_orbit(2); periodic_orbit(3); periodic_orbit(4)];
for i=1:k    
    tn = periodic_orbit(1+4*(i-1));
    xn = periodic_orbit(2+4*(i-1));
    vn = periodic_orbit(3+4*(i-1));
    wn = periodic_orbit(4+4*(i-1));
    tn1 = periodic_orbit(1+4*i);    
% The reset law's output
     R_out(1+2*(i-1):2+2*(i-1)) = R_map(tn, xn, vn, wn, a, E, e);
% the reset law's jacobian, dv = [dvdt, dvdx, dvdv, dvdw], dw = [dwdt, dwdx, dwdv, dwdw] 
    dv = [(8*E*a*pi^2*xn*sin(2*pi*tn)*(e+1))/(1+(2*a*xn)^2), ...
          (-2*a*(e+1)*(4*a^2*wn*xn^2-wn+2*pi*E*cos(2*pi*tn)+4*a*vn*xn-8*E*a^2*xn^2*pi*cos(2*pi*tn)))/((1+(2*a*xn)^2)^2), ...
          (1-4*e*a^2*xn^2)/(1+(2*a*xn)^2), 2*a*xn*(e+1)/(1+(2*a*xn)^2)];
    dw = [(-4*E*pi^2*sin(2*pi*tn)*(e+1))/(1+(2*a*xn)^2), (2*a*(e+1)*(vn-4*a^2*vn*xn^2+4*a*wn*xn-8*E*a*xn*pi*cos(2*pi*tn)))/((1+(2*a*xn)^2)^2), ...
           (2*a*xn*(e+1))/(1+(2*a*xn)^2), (4*a^2*xn^2-e)/(1+(2*a*xn)^2)];      
% Constructing the zero-finding problem's Jacobian
% partial of tn
    block1(1+4*(i-1),1+4*(i-1)) = 2*a*(xn+R_out(1+2*(i-1))*(tn1-tn))*(dv(1)*(tn1-tn)-R_out(1+2*(i-1))) - 2*pi*E*cos(2*pi*tn)+R_out(2+2*(i-1))-(dw(1)+9.81)*(tn1-tn);
    block1(2+4*(i-1),1+4*(i-1)) = dv(1)*(tn1-tn)-R_out(1+2*(i-1));
    block1(3+4*(i-1),1+4*(i-1)) = dv(1);    
    block1(4+4*(i-1),1+4*(i-1)) = dw(1)+9.81;    
% partial of xn
    block1(1+4*(i-1),2+4*(i-1)) = 2*a*(tn1-tn)*(R_out(1+2*(i-1))+xn*dv(2))+2*a*dv(2)*R_out(1+2*(i-1))*(tn1-tn)^2-dw(2)*(tn1-tn);
    block1(2+4*(i-1),2+4*(i-1)) = 1+dv(2)*(tn1-tn);
    block1(3+4*(i-1),2+4*(i-1)) = dv(2);
    block1(4+4*(i-1),2+4*(i-1)) = dw(2);     
% partial of vn
    block1(1+4*(i-1),3+4*(i-1)) = 2*a*dv(3)*(tn1-tn)*(xn+R_out(1+2*(i-1))*(tn1-tn))-dw(3)*(tn1-tn);
    block1(2+4*(i-1),3+4*(i-1)) = dv(3)*(tn1-tn);    
    block1(3+4*(i-1),3+4*(i-1)) = dv(3);
    block1(4+4*(i-1),3+4*(i-1)) = dw(3);    
% partial of wn
    block1(1+4*(i-1),4+4*(i-1)) = 2*a*dv(4)*(tn1-tn)*(xn+R_out(1+2*(i-1))*(tn1-tn))-dw(4)*(tn1-tn);
    block1(2+4*(i-1),4+4*(i-1)) = dv(4)*(tn1-tn);    
    block1(3+4*(i-1),4+4*(i-1)) = dv(4);
    block1(4+4*(i-1),4+4*(i-1)) = dw(4);    
% partials of tn+1
    block2(1+4*(i-1),1+4*(i-1)) = 2*a*R_out(1+2*(i-1))*(xn+R_out(1+2*(i-1))*(tn1-tn))+2*pi*E*cos(2*pi*tn1)-R_out(2+2*(i-1))+9.81*(tn1-tn);
    block2(2+4*(i-1),1+4*(i-1)) = R_out(1+2*(i-1));   
    block2(3+4*(i-1),1+4*(i-1)) = 0;
    block2(4+4*(i-1),1+4*(i-1)) = -9.81;    
% partials of xn+1
    block2(1+4*(i-1),2+4*(i-1)) = 0;
    block2(2+4*(i-1),2+4*(i-1)) = -1;    
    block2(3+4*(i-1),2+4*(i-1)) = 0;
    block2(4+4*(i-1),2+4*(i-1)) = 0;    
% partials of vn+1
    block2(1+4*(i-1),3+4*(i-1)) = 0;
    block2(2+4*(i-1),3+4*(i-1)) = 0;    
    block2(3+4*(i-1),3+4*(i-1)) = -1;
    block2(4+4*(i-1),3+4*(i-1)) = 0;    
% partials of wn+1
    block2(1+4*(i-1),4+4*(i-1)) = 0;
    block2(2+4*(i-1),4+4*(i-1)) = 0;  
    block2(3+4*(i-1),4+4*(i-1)) = 0;    
    block2(4+4*(i-1),4+4*(i-1)) = -1;
end
y = block1;
for j=1:k-1
    y(1+4*(j-1):4+4*(j-1), 1+4*j:4+4*j) = block2(1+4*(j-1):4+4*(j-1),1+4*(j-1):4+4*(j-1)); 
end
y(1+4*(k-1):4+4*(k-1), 1:4) = block2(1+4*(k-1):4+4*(k-1),1+4*(k-1):4+4*(k-1));
end