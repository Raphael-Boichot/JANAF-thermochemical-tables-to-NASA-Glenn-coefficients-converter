function [f2] = dist_enthalpie(x,t,y_exp)
global a1 a2 a3 a4 a5;
a6=x(1) ;

m = length(t) ;
p = zeros(m,1) ;
for k = 1:m
   tk = t(k) ;
   p(k) = a1+a2.*tk./2+a3.*(tk.^2./3)+a4.*(tk.^3./4)+a5.*(tk.^4./5)+a6./tk ;
end

r = p - y_exp ;
f2 = sum(r.^2) ;

