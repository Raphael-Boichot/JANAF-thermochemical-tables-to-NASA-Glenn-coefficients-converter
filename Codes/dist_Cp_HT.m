function [f2] = dist_Cp_HT(x,t,y_exp)
global Cp_R_Cutting;
global Cp_R_slope_Cutting;
global Cutting_temperature;
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

r = p - y_exp ;%residuals
Cp_cutting=a1+a2.*Cutting_temperature+a3.*Cutting_temperature.^2+a4.*Cutting_temperature.^3+a5.*Cutting_temperature.^4;
Cp_cutting_plus=a1+a2.*(Cutting_temperature+1)+a3.*(Cutting_temperature+1).^2+a4.*(Cutting_temperature+1).^3+a5.*(Cutting_temperature+1).^4;
Cp_R_slope_Cutting_local=Cp_cutting_plus-Cp_cutting;%here we get the local Cp and Cp derivative at BT/HT transition
%here we enforce Cp slope and Cp value continuity
f2 = sum(r.^2)+10*abs(Cp_R_Cutting-p(1))+10*abs(Cp_R_slope_Cutting_local-Cp_R_slope_Cutting);
