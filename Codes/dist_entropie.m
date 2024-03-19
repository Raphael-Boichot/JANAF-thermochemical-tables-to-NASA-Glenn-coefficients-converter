function [f2] = dist_entropie(x,t,y_exp)
global a1 a2 a3 a4 a5;
a7=x(1) ;

m = length(t) ;
p = zeros(m,1) ;
for k = 1:m
   tk = t(k) ;
   p(k) = a1.*log(tk)+a2.*tk+a3.*(tk.^2./2)+a4.*(tk.^3./3)+a5.*(tk.^4./4)+a7 ;
end

r = p - y_exp ;
f2 = sum(r.^2) ;

