function [out] = S_map(t, s, x, v, w, a, E)
% This function is the S map, it takes a 5D vector [tn+1, tn, xn, R(tn, xn, vn, wn)] 
% and outputs a 4D vector [xn+1, yn+1, vn+1, wn+1].
t = t(:).'; % columnn vector
[out] = [x + v.*(t - s); 
         a.*x.^2 + E.*sin(2.*pi.*s) + w.*(t - s) - 4.905*(t - s).^2; 
         v.*(1+t-t); 
         w - 9.81.*(t - s)];
% The difference in the third entry "t-t" is to make sure the function
% outputs a valid size matrix when the input "t" is a vector.
end
