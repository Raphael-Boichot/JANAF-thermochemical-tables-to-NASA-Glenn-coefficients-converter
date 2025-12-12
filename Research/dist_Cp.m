function [f2] = dist_Cp(x,t,y_exp)
a1 = x(1) ;
a2 = x(2) ;
a3 = x(3) ;
a4 = x(4) ;
a5 = x(5) ;

m = length(t) ;
p = zeros(m,1) ;
for k = 1:m
   tk = t(k) ;
   p(k) = a1+a2.*tk+a3.*tk.^2+a4.*tk.^3+a5.*tk.^4 ;
end

r = p - y_exp ;
f2 = sum(r.^2);

