%Created by Raphael BOICHOT around 2018, updated in March 2024 for public release
%Format of output data : T(K) H/(RT)(-) Cp/R(-) S/R(-) in columns
clc;
clear;
close all;
R=8.314462618;%(J/kmol.K)
cal=4.184;
R_cal=R/cal;
Tref = 298.15;
Default_fontsize=16;
fileID = fopen('Input.dat','r');
A = fread(fileID,45,'*char');
Data=cell2mat(textscan(char(fread(fileID,35,'*char')),'%f'));
Low_temp = Data(1);
High_temp = Data(2);
Cutting_temp = round(Data(3));
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
T_BT=[Tref,Low_temp:1:Cutting_temp];
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
% raw_data=[T H_RT Cp_R S_R ];
% save thermo_data.txt raw_data -ascii;

fprintf('Cutting_temp = %g\n', Cutting_temp);
fprintf('Low coeffs a1..a7 =\n');
fprintf('%g ', a1_BT,a2_BT,a3_BT,a4_BT,a5_BT,a6_BT,a7_BT);
fprintf('\nHigh coeffs a1..a7 =\n');
fprintf('%g ', a1_HT,a2_HT,a3_HT,a4_HT,a5_HT,a6_HT,a7_HT);
fprintf('\n');

% compute H/RT at 298.15 using both sets
deltaHf = R*Tref * ( a1_BT + a2_BT*Tref/2 + a3_BT*Tref^2/3 + ...
               a4_BT*Tref^3/4 + a5_BT*Tref^4/5 + a6_BT/Tref );
fprintf('H(298.15) using LOW coefficients  = %g J/mol\n', deltaHf);

%hypothèses prises pour Co(acac)2
%On part des poynomes NASA du ligand acac seul
%La molécule complète a un Cp et un S de 2x Cp et S du ligand acac %(https://rmg.mit.edu/database/thermo/libraries/CHO/106/)
%Le deltaH de formation de acac et Co(acac)3 sont connus (Inorg.Chem.2025,64,21405−21418)
%On intégre H à partir du Cp(T) avec le deltaH de formation réel de Co(acac)2 (2/3) de celui de Co(acac)3

deltaHf=-267.6; %kcal.mol-1, valeur réelle pour Co(acac3)
Nb_ligand=2;
deltaHf_new=deltaHf*1000*cal*(2/3) %J/moles, 2/3 pour Co(acac)2
S_R_new=Nb_ligand*S_R;
Cp_R_new=Nb_ligand*Cp_R;

%recréation de H/RT from scratch à partir de Cp et de deltaH
Cp_fun_BT = @(T_int) Nb_ligand*R*(a1_BT + a2_BT.*T_int + a3_BT.*T_int.^2 + a4_BT.*T_int.^3 + a5_BT.*T_int.^4);
Cp_fun_HT = @(T_int) Nb_ligand*R*(a1_HT + a2_HT.*T_int + a3_HT.*T_int.^2 + a4_HT.*T_int.^3 + a5_HT.*T_int.^4);
cut=find(T==Cutting_temp)
for index=1:1:cut(1)
T_int=T(index);
H_new(index,1)=deltaHf_new+integral(Cp_fun_BT,298.15,T_int);
end
for index=cut(2):1:length(T)
T_int=T(index);
H_new(index,1)=deltaHf_new+integral(Cp_fun_HT,298.15,T_int);
end
%force la continuité de H
offset=H_new(cut(1))-H_new(cut(2))
H_new(cut(2):end)=H_new(cut(2):end)+offset;
H_RT_new=H_new./(R*T);

figure(4)
plot(T,H_RT_new.*T*R,'LineWidth',2)
title('New enthalpy vs T','Fontsize',Default_fontsize)
ylabel('H (J/mol)','Fontsize',Default_fontsize)
xlabel('Temperature in K','Fontsize',Default_fontsize);
set(gca,'FontSize',Default_fontsize)
saveas(gcf,'H_R_NASA_new.png');

figure(5)
plot(T,Cp_R_new*(R),'LineWidth',2)
title('New heat capacity vs T','Fontsize',Default_fontsize)
ylabel('Cp (J/(mol.K))','Fontsize',Default_fontsize)
xlabel('Temperature in K','Fontsize',Default_fontsize);
set(gca,'FontSize',Default_fontsize)
saveas(gcf,'Cp_R_NASA_new.png');

figure(6)
plot(T,S_R_new*(R),'LineWidth',2)
title('New entropy vs T','Fontsize',Default_fontsize)
ylabel('S (J/(mol.K))','Fontsize',Default_fontsize)
xlabel('Temperature in K','Fontsize',Default_fontsize);
set(gca,'FontSize',Default_fontsize)
saveas(gcf,'S_R_NASA_new.png');


[T_unique, idx_unique] = unique(T);
T = T(idx_unique); 
H_RT_new = H_RT_new(idx_unique);          
Cp_R_new = Cp_R_new(idx_unique);   
S_R_new = S_R_new(idx_unique);   

raw_data=[T H_RT_new Cp_R_new S_R_new];
save JANAF_input.txt raw_data -ascii;