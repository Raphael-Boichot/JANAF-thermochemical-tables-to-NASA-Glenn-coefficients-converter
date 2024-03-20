%Created by Raphael BOICHOT around 2018, updated in March 2024 for public release
%Format of entry data : T(K) H/(RT)(-) Cp/R(-) S/R(-) in columns
clc;
clear;


Low_temp_threshold=800;    %Low temperature threshold for cutting temperature search
High_temp_threshold=1200;  %High temperature threshold for cutting temperature search
Searching_T_step=10;       %Searching step size in K for cutting temperature
n=5;                       %Skip data every n points to accelerate algorithm

R=8.314;
Default_font_size=14;
format long
global a1 a2 a3 a4 a5 Cp_R_Cutting Cp_R_slope_Cutting Cutting_temperature;
raw_data = load('JANAF_input.txt') ;
T_exp = raw_data(:,1);
H_R_exp = raw_data(:,2);
Cp_R_exp = raw_data(:,3);
S_R_exp = raw_data(:,4);

OPTIONS = optimset(@fminsearch) ;
OPTIONS = optimset(OPTIONS,'MaxFunEvals',4000,'MaxIter',4000,'TolFun',1e-6,'TolX',1e-6) ;
OPTIONS = optimset(OPTIONS,'Display','off') ;
k=0;
disp('Searching the best Low temp/High temp cutting point for polynomials')
tic
for Cutting_temperature=Low_temp_threshold:Searching_T_step:High_temp_threshold
    k=k+1;
    disp(['Progression: ',num2str((Cutting_temperature-Low_temp_threshold)/(High_temp_threshold-Low_temp_threshold)*100,3),'%'])
    cutting_pos=find(T_exp==Cutting_temperature);
    %******************************************************************************
    %We began by low temperature range: optimization without constraints
    T_BT_exp=T_exp(1:n:cutting_pos);
    H_R_BT_exp=H_R_exp(1:n:cutting_pos);
    S_R_BT_exp=S_R_exp(1:n:cutting_pos);
    Cp_R_BT_exp=Cp_R_exp(1:n:cutting_pos);
    x_0 = [10 1e-4 1e-7 1e-11 1e-13];%some realistic starting values
    %first the polynomials for Cp (a1:a5)
    [x_opt, res0] = fminsearch(@dist_Cp,x_0,OPTIONS,T_BT_exp,Cp_R_BT_exp) ;
    %We need these values for the H and S optimization
    a1=x_opt(1);
    a2=x_opt(2);
    a3=x_opt(3);
    a4=x_opt(4);
    a5=x_opt(5);
    %Then a6 for H, a1:a5 been fixed now (H is the Cp intergral)
    x_0 = 1;%some random starting value
    [x_opt,res1] = fminsearch(@dist_enthalpie,x_0,OPTIONS,T_BT_exp,H_R_BT_exp) ;
    a6_BT=x_opt;
    x_0 = 1;%some random starting value
    %Then a7 for S, a1:a6 been fixed now (S is the H intergral over T)
    [x_opt,res2] = fminsearch(@dist_entropie,x_0,OPTIONS,T_BT_exp,S_R_BT_exp) ;
    a7_BT=x_opt;

    %Then we deal with the high temperature range: optimization with
    %constraints of Cp slope and Cp continuity
    T_HT_exp=T_exp(cutting_pos:n:end);
    H_R_HT_exp=H_R_exp(cutting_pos:n:end);
    S_R_HT_exp=S_R_exp(cutting_pos:n:end);
    Cp_R_HT_exp=Cp_R_exp(cutting_pos:n:end);
    x_0 = [10 1e-4 1e-7 1e-11 1e-13];%some realistic starting values
    %here we get the local Cp and Cp derivative at BT/HT transition
    T=Cutting_temperature;
    Cp_R_Cutting=a1+a2*T+a3*T.^2+a4.*T.^3+a5*T.^4;
    T=T+1;%local slope
    Cp_R_Cutting_plus=a1+a2*T+a3*T.^2+a4*T.^3+a5*T.^4;
    Cp_R_slope_Cutting=Cp_R_Cutting_plus-Cp_R_Cutting;

    %first the polynomials for Cp (a1:a5) with forced Cp value and slope at
    %the HT/BT boundary
    [x_opt,res3] = fminsearch(@dist_Cp_HT,x_0,OPTIONS,T_HT_exp,Cp_R_HT_exp) ;
    %We need these values for the H and S optimization
    a1=x_opt(1);
    a2=x_opt(2);
    a3=x_opt(3);
    a4=x_opt(4);
    a5=x_opt(5);
    x_0 = 1 ;%some random starting value
    [x_opt,res4] = fminsearch(@dist_enthalpie,x_0,OPTIONS,T_HT_exp,H_R_HT_exp) ;
    a6_HT=x_opt;
    x_0 = 1 ;%some random starting value
    [x_opt,res5] = fminsearch(@dist_entropie,x_0,OPTIONS,T_HT_exp,S_R_HT_exp) ;
    a7_HT=x_opt;
    %******************************************************************************
    res(k)=res0+res1+res2+res3+res4+res5;
end
toc
[M,I]=min(res);
A=Low_temp_threshold:Searching_T_step:High_temp_threshold;
disp(['Best cutting temperature = ',num2str(A(I))]);
disp(['Residuals = ',num2str(res(I))]);
disp(' ');

%********************Polynomial optimization********************************
Low_temp_threshold=T_exp(1);
Cutting_temperature=A(I);
High_temp_threshold=T_exp(end);
cutting_pos=find(T_exp==Cutting_temperature);

%We began by low temperature range: optimization without constraints
T_BT_exp=T_exp(1:n:cutting_pos);
H_R_BT_exp=H_R_exp(1:n:cutting_pos);
S_R_BT_exp=S_R_exp(1:n:cutting_pos);
Cp_R_BT_exp=Cp_R_exp(1:n:cutting_pos);
x_0 = [10 1e-4 1e-7 1e-11 1e-13];%some realistic starting values
%first the polynomials for Cp (a1:a5)
[x_opt, res0] = fminsearch(@dist_Cp,x_0,OPTIONS,T_BT_exp,Cp_R_BT_exp) ;
a1_BT=x_opt(1);
a2_BT=x_opt(2);
a3_BT=x_opt(3);
a4_BT=x_opt(4);
a5_BT=x_opt(5);
%We need these values for the H and S optimization
a1=x_opt(1);
a2=x_opt(2);
a3=x_opt(3);
a4=x_opt(4);
a5=x_opt(5);
%Then a6 for H, a1:a5 been fixed now (H is the Cp intergral)
x_0 = 1;%some random starting value
[x_opt,res1] = fminsearch(@dist_enthalpie,x_0,OPTIONS,T_BT_exp,H_R_BT_exp) ;
a6_BT=x_opt;
x_0 = 1;%some random starting value
%Then a7 for S, a1:a6 been fixed now (S is the H intergral over T)
[x_opt,res2] = fminsearch(@dist_entropie,x_0,OPTIONS,T_BT_exp,S_R_BT_exp) ;
a7_BT=x_opt;
%Then we deal with the high temperature range: optimization with
%constraints of Cp slope and Cp continuity
T_HT_exp=T_exp(cutting_pos:n:end);
H_R_HT_exp=H_R_exp(cutting_pos:n:end);
S_R_HT_exp=S_R_exp(cutting_pos:n:end);
Cp_R_HT_exp=Cp_R_exp(cutting_pos:n:end);
x_0 = [10 1e-4 1e-7 1e-11 1e-13];%some realistic starting values

%here we get the local Cp and Cp derivative at BT/HT transition
T=Cutting_temperature;
Cp_R_Cutting=a1_BT+a2_BT.*T+a3_BT.*T.^2+a4_BT.*T.^3+a5_BT.*T.^4;
T=T+1;%local slope
Cp_R_Cutting_plus=a1_BT+a2_BT.*T+a3_BT.*T.^2+a4_BT.*T.^3+a5_BT.*T.^4;
Cp_R_slope_Cutting=Cp_R_Cutting_plus-Cp_R_Cutting;

[x_opt,res3] = fminsearch(@dist_Cp_HT,x_0,OPTIONS,T_HT_exp,Cp_R_HT_exp) ;
a1_HT=x_opt(1);
a2_HT=x_opt(2);
a3_HT=x_opt(3);
a4_HT=x_opt(4);
a5_HT=x_opt(5);
%We need these values for the H and S optimization
a1=x_opt(1);
a2=x_opt(2);
a3=x_opt(3);
a4=x_opt(4);
a5=x_opt(5);
x_0 = 1 ;%some random starting value
[x_opt,res4] = fminsearch(@dist_enthalpie,x_0,OPTIONS,T_HT_exp,H_R_HT_exp) ;
a6_HT=x_opt;
x_0 = 1 ;
[x_opt,res5] = fminsearch(@dist_entropie,x_0,OPTIONS,T_HT_exp,S_R_HT_exp) ;
a7_HT=x_opt;%some random starting value

%data for plotting
T=T_BT_exp;
Cp_R_BT=a1_BT+a2_BT.*T+a3_BT.*T.^2+a4_BT.*T.^3+a5_BT.*T.^4;
H_R_BT=a1_BT+a2_BT.*T/2+a3_BT.*T.^2/3+a4_BT.*T.^3/4+a5_BT.*T.^4/5+a6_BT./T;
S_R_BT=a1_BT.*log(T)+a2_BT.*T+a3_BT.*T.^2/2+a4_BT.*T.^3/3+a5_BT.*T.^4/4+a7_BT;

%data for plotting
T=T_HT_exp;
Cp_R_HT=a1_HT+a2_HT.*T+a3_HT.*T.^2+a4_HT.*T.^3+a5_HT.*T.^4;
H_R_HT=a1_HT+a2_HT.*T/2+a3_HT.*T.^2/3+a4_HT.*T.^3/4+a5_HT.*T.^4/5+a6_HT./T;
S_R_HT=a1_HT.*log(T)+a2_HT.*T+a3_HT.*T.^2/2+a4_HT.*T.^3/3+a5_HT.*T.^4/4+a7_HT;

figure('Position',[100 100 1200 1000]);
figure(1)
subplot(2,2,1)
plot(A,log10(res),A(I),log10(res(I)),'o','MarkerFaceColor','red','LineWidth',2);
title('Residuals vs T cutting','Fontsize',Default_font_size)
ylabel('Log10 residuals','Fontsize',Default_font_size)
xlabel('Temperature(K)','Fontsize',Default_font_size);
text(A(I)+25,log10(res(I)),'Local minimum');
set(gca,'FontSize',Default_font_size)

subplot(2,2,2)
plot(T_BT_exp,H_R_BT.*T_BT_exp*(R),'r-',T_BT_exp,H_R_BT_exp.*T_BT_exp.*R,'--g',...
    T_HT_exp,H_R_HT.*T_HT_exp*(R),'r-',T_HT_exp,H_R_HT_exp.*T_HT_exp.*R,'--g','LineWidth',2)
title('Enthalpy vs T','Fontsize',Default_font_size)
ylabel('H (J/mol)','Fontsize',Default_font_size)
xlabel('Temperature(K)','Fontsize',Default_font_size);
legend('Fitted data','Input data','Location','southeast');
set(gca,'FontSize',Default_font_size)

subplot(2,2,3)
plot(T_BT_exp,Cp_R_BT*(R),'r-',T_BT_exp,Cp_R_BT_exp.*R,'--g',...
    T_HT_exp,Cp_R_HT*(R),'r-',T_HT_exp,Cp_R_HT_exp.*R,'--g','LineWidth',2);
title('Heat capacity vs T','Fontsize',Default_font_size)
ylabel('Cp (J/(mol.K))','Fontsize',Default_font_size)
xlabel('Temperature(K)','Fontsize',Default_font_size);
legend('Fitted data','Input data','Location','southeast');
set(gca,'FontSize',Default_font_size)

subplot(2,2,4)
plot(T_BT_exp,S_R_BT*(R),'r-',T_BT_exp,S_R_BT_exp.*R,'--g',...
    T_HT_exp,S_R_HT*(R),'r-',T_HT_exp,S_R_HT_exp.*R,'--g','LineWidth',2);
title('Entropy vs T','Fontsize',Default_font_size)
ylabel('S (J/(mol.K))','Fontsize',Default_font_size)
xlabel('Temperature in K','Fontsize',Default_font_size);
legend('Fitted data','Input data','Location','southeast');
set(gca,'FontSize',Default_font_size)

saveas(gcf,'NASA_fitting.png');

f='%+10.8e\n';
%****************Free comments, must fit within columns 1:24***
column_1_24='CO2               L 7/88';
%*************A***NA***NA***NA***N atom formula: Atom***Number, and so on...
%************must fit within columns 25:44
column_25_44='C   1O   2    0    0';
% G for gaseous, C for condensed, must be at column 45
column_45='G';
% All the next text must fit EXACTLY within the same columns/lines
disp('****************start of text to copy paste*************************************');
disp([column_1_24,column_25_44,column_45,'    ',...
    num2str(Low_temp_threshold,'%.2f'),' ',num2str(High_temp_threshold,'%.2f'),' ',num2str(Cutting_temperature,'%.2f'),'        1']);
disp([num2str(a1_HT,f),num2str(a2_HT,f),num2str(a3_HT,f),num2str(a4_HT,f),num2str(a5_HT,f),'    2']);
disp([num2str(a6_HT,f),num2str(a7_HT,f),num2str(a1_BT,f),num2str(a2_BT,f),num2str(a3_BT,f),'    3']);
disp([num2str(a4_BT,f),num2str(a5_BT,f),num2str(a6_BT,f),num2str(a7_BT,f),'                   4']);
disp('*****************end of text to copy paste**************************************');

fileID = fopen('Output.dat','w');
fwrite(fileID,[column_1_24,column_25_44,column_45,'    ',...
    num2str(Low_temp_threshold,'%.2f'),'   ',num2str(High_temp_threshold,'%.2f'),' ',num2str(Cutting_temperature,'%.2f'),'    ','0 1']);
fwrite(fileID,char(13));
fwrite(fileID,newline);
fwrite(fileID,[num2str(a1_HT,'%+10.8e\n'),num2str(a2_HT,'%+10.8e\n'),num2str(a3_HT,'%+10.8e\n'),num2str(a4_HT,'%+10.8e\n'),num2str(a5_HT,'%+10.8e\n'),'    2']);
fwrite(fileID,char(13));
fwrite(fileID,newline);
fwrite(fileID,[num2str(a6_HT,'%+10.8e\n'),num2str(a7_HT,'%+10.8e\n'),num2str(a1_BT,'%+10.8e\n'),num2str(a2_BT,'%+10.8e\n'),num2str(a3_BT,'%+10.8e\n'),'    3']);
fwrite(fileID,char(13));
fwrite(fileID,newline);
fwrite(fileID,[num2str(a4_BT,'%+10.8e\n'),num2str(a5_BT,'%+10.8e\n'),num2str(a6_BT,'%+10.8e\n'),num2str(a7_BT,'%+10.8e\n'),'                   4']);
fwrite(fileID,char(13));
fwrite(fileID,newline);
fclose(fileID);
