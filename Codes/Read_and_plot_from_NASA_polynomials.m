%Created by Raphael BOICHOT around 2018, updated in March 2024 for public release
%Format of output data : T(K) H/(RT)(-) Cp/R(-) S/R(-) in columns
clc;
clear;
R=8.314462618;%(J/kmol.K)
Default_fontsize=16;
fileID = fopen('Input.dat','r');
A = fread(fileID,45,'*char');
Data=cell2mat(textscan(char(fread(fileID,35,'*char')),'%f'));
Low_temp = Data(1);
High_temp = Data(2);
Cutting_temp = Data(3);
fread(fileID,2,'*char');
a1_HT = cell2mat(textscan(char(fread(fileID,15,'*char')),'%f'));
a2_HT = cell2mat(textscan(char(fread(fileID,15,'*char')),'%f'));
a3_HT = cell2mat(textscan(char(fread(fileID,15,'*char')),'%f'));
a4_HT = cell2mat(textscan(char(fread(fileID,15,'*char')),'%f'));
a5_HT = cell2mat(textscan(char(fread(fileID,15,'*char')),'%f'));
fread(fileID,7,'*char');
a6_HT = cell2mat(textscan(char(fread(fileID,15,'*char')),'%f'));
a7_HT = cell2mat(textscan(char(fread(fileID,15,'*char')),'%f'));
a1_BT = cell2mat(textscan(char(fread(fileID,15,'*char')),'%f'));
a2_BT = cell2mat(textscan(char(fread(fileID,15,'*char')),'%f'));
a3_BT = cell2mat(textscan(char(fread(fileID,15,'*char')),'%f'));
fread(fileID,7,'*char');
a4_BT = cell2mat(textscan(char(fread(fileID,15,'*char')),'%f'));
a5_BT = cell2mat(textscan(char(fread(fileID,15,'*char')),'%f'));
a6_BT = cell2mat(textscan(char(fread(fileID,15,'*char')),'%f'));
a7_BT = cell2mat(textscan(char(fread(fileID,15,'*char')),'%f'));
fclose(fileID);

%Plotting
T_BT=Low_temp:1:Cutting_temp;
H_R_BT=a1_BT+a2_BT.*T_BT/2+a3_BT.*T_BT.^2/3+a4_BT.*T_BT.^3/4+a5_BT.*T_BT.^4/5+a6_BT./T_BT;
Cp_R_BT=a1_BT+a2_BT.*T_BT+a3_BT.*T_BT.^2+a4_BT.*T_BT.^3+a5_BT.*T_BT.^4;
S_R_BT=a1_BT.*log(T_BT)+a2_BT.*T_BT+a3_BT.*T_BT.^2/2+a4_BT.*T_BT.^3/3+a5_BT.*T_BT.^4/4+a7_BT;

T_HT=Cutting_temp:1:High_temp;
H_R_HT=a1_HT+a2_HT.*T_HT/2+a3_HT.*T_HT.^2/3+a4_HT.*T_HT.^3/4+a5_HT.*T_HT.^4/5+a6_HT./T_HT;
Cp_R_HT=a1_HT+a2_HT.*T_HT+a3_HT.*T_HT.^2+a4_HT.*T_HT.^3+a5_HT.*T_HT.^4;
S_R_HT=a1_HT.*log(T_HT)+a2_HT.*T_HT+a3_HT.*T_HT.^2/2+a4_HT.*T_HT.^3/3+a5_HT.*T_HT.^4/4+a7_HT;

figure(1)
plot(T_BT,H_R_BT.*T_BT*(R),T_HT,H_R_HT.*T_HT*R,'LineWidth',2)
title('Enthalpy vs T','Fontsize',Default_fontsize)
ylabel('H (J/mol)','Fontsize',Default_fontsize)
xlabel('Temperature in K','Fontsize',Default_fontsize);
set(gca,'FontSize',Default_fontsize)
saveas(gcf,'H_R_NASA.png');

figure(2)
plot(T_BT,Cp_R_BT*(R),T_HT,Cp_R_HT*R,'LineWidth',2)
title('Heat capacity vs T','Fontsize',Default_fontsize)
ylabel('Cp (J/(mol.K))','Fontsize',Default_fontsize)
xlabel('Temperature in K','Fontsize',Default_fontsize);
set(gca,'FontSize',Default_fontsize)
saveas(gcf,'Cp_R_NASA.png');

figure(3)
plot(T_BT,S_R_BT*(R),T_HT,S_R_HT*R,'LineWidth',2)
title('Entropy vs T','Fontsize',Default_fontsize)
ylabel('S (J/(mol.K))','Fontsize',Default_fontsize)
xlabel('Temperature in K','Fontsize',Default_fontsize);
set(gca,'FontSize',Default_fontsize)
saveas(gcf,'S_R_NASA.png');

T=[T_BT T_HT]';
S_R=[S_R_BT S_R_HT]';
H_RT=[H_R_BT H_R_HT]';
Cp_R=[Cp_R_BT Cp_R_HT]';
raw_data=[T H_RT Cp_R S_R ];
save thermo_data.txt raw_data -ascii;

Tref = 298.15;

fprintf('Cutting_temp = %g\n', Cutting_temp);
fprintf('Low coeffs a1..a7 =\n');
fprintf('%g ', a1_BT,a2_BT,a3_BT,a4_BT,a5_BT,a6_BT,a7_BT);
fprintf('\nHigh coeffs a1..a7 =\n');
fprintf('%g ', a1_HT,a2_HT,a3_HT,a4_HT,a5_HT,a6_HT,a7_HT);
fprintf('\n');

% compute H/RT at 298.15 using both sets
H_R_low_ref = a1_BT + a2_BT*Tref/2 + a3_BT*Tref^2/3 + a4_BT*Tref^3/4 + a5_BT*Tref^4/5 + a6_BT/Tref;
H_R_high_ref= a1_HT + a2_HT*Tref/2 + a3_HT*Tref^2/3 + a4_HT*Tref^3/4 + a5_HT*Tref^4/5 + a6_HT/Tref;

H_low = H_R_low_ref * R * Tref;
H_high= H_R_high_ref * R * Tref;

fprintf('H(298) using LOW coefficients  = %g J/mol\n', H_low);
fprintf('H(298) using HIGH coefficients = %g J/mol\n', H_high);