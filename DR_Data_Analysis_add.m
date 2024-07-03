function [Alpha0_All,Alpha_All,Cell_DR_Counts_DCCT_Data,filename]=DR_Data_Analysis_add(Z_num,Am,qcharge,Eion,Ie,deltan,nummax,deltamn1,cyclenum1,cycletime1,type,pick)
% [Alpha_All,filename,Cell_DR_Counts_DCCT_Data]=DR_Data_Analysis(deltan,deltamn,nummax,cyclenum,cycletime)
% This function writen aims for analysis on the DR experiment data, obtain
% the first spectra of RR and DR, and help us judge the operation of DR
% experiment is right. The meanings of parameters are in the following.

% Z_num: the atomic number of ion

% Am: the mass number of ion

% Eion: the kinetic energy per nuclear of ion (unit: Mev/u)

% qcharge: the charge number of ion

% Ie: the current of electron beam

% deltan: the detuning counts data number you want to calculate, depends on 
% the detuning time you taken. For example, for detuning time of 10/20 ms, 
% recommended number of deltan maybe be 6-8-10/15-18-20.

% nummax: the maximum number we judge the location of NIM2 signal is
% right, by which we judge what range we pick up detuning counts data. For 
% example, for the detuning time of 10/20 ms, recommended number of nummax 
% Should be large than 10[such as 15]/20[such as 25].

% deltamn: the cooling counts data number you want to calculate, depends on 
% the cooling time you taken. For example, for cooling time of 90/80 ms,
% recommended number of deltamn maybe be 80-90/70-80.

% cyclenum: the sample time cycle number, i.e., the sum time number, 
% which equals to the detuning time number adding to the cooling time
% number. For example a 100 ms sample time cycle, which equals to the
% detuning time number of 10 adds to the cooling time number of 90.

% cycletime: the search DCCT counts cycle time number [99, recommended number 90]

% type:if equcals to 1 ,do without calculating  Alpha_d/Alpha_c and add normalization.

% pick:decide how many points of detuning counts you want to skip. For
% example,if you wants to start calculate at the second point of detuning,pick should
% be 1.

% Author: WQ Xu    $ Company: USTC  

% clear desktop
% tic;  % Obtain the elapsed time of this program

clc;
close all;

%% Open DR Data files directory and select data files
[filename,pathname,flag]=uigetfile({'*','All Files'},...
    'Please Slecte DR Data Files','MultiSelect','on');
if ~flag
    return;
end

filename=cellstr(filename);
nfiles=length(filename);
filename=filename';

% Clear the unuseful data filenames in the worksapce of MATLAB, such as: 
% a) two repeated files with the same detuning voltage and detuning pulse 
% width due to some reasons, such as mistaken operations. b) for more long
% time accumulation counts, etc. By the time series, we maybe delete the 
% former file. 
% Other situations maybe handle depending on the situations.

%% Parameters of electron and ion
me=9.1e-31;mion=Am*1.67e-27+me*(Z_num-qcharge);           % in unit of Kg
c=3e8;e=1.6e-19;                           % in unit of SCI
lc=3.4;

% Obtain the detuning devoltages correspanding to the data files
detuning_voltages=zeros(nfiles,2,'double');

for a=1:nfiles
    
    num=length(filename{a});
    for b=1:num
      if strcmp(filename{a}(b),'V')
         detuning_voltages(a,1)=str2num(filename{a}(1:b-1));
         detuning_voltages(a,2)=-detuning_voltages(a,1);
      end
    end
end


epsinon=8.854187818e-12;
% Initial kinetic energy per nuclear of ion [Eion]
gamaion=1+Eion/931.5;
betaion=sqrt(1-1/gamaion^2);
rad_tube=25/2;
radius_e=2.95;

% Calculation of space-charge effect for electron beam
%Usp0=(gamaion-1-Ucath/0.511/1e3)*0.511*1e6;
Usp0=-0.93*(Ie/(4*pi*epsinon*betaion*c))*(1+2*log(rad_tube/radius_e));
Para=-Usp0*betaion;
Ucath=(gamaion-1-Usp0/511/1e3)*0.511*1e6;

%=================Modify the relative energy================
gamaee=zeros(nfiles,2,'double');
gamaee(:,1)=gamaion-(1e-6)*detuning_voltages(:,1)/0.511;
gamaee(:,2)=gamaion-(1e-6)*detuning_voltages(:,2)/0.511;

betaee=zeros(nfiles,2,'double');
betaee(:,1)=sqrt(1-1./gamaee(:,1).^2);
betaee(:,2)=sqrt(1-1./gamaee(:,2).^2);

Usp=zeros(nfiles,2,'double');
Usp(:,1)=-Para./betaee(:,1);
Usp(:,2)=-Para./betaee(:,2);

gamae=zeros(nfiles,2,'double');
gamae(:,1)=1+(1e-6)*(Ucath-detuning_voltages(:,1)+Usp(:,1))/0.511;
gamae(:,2)=1+(1e-6)*(Ucath-detuning_voltages(:,2)+Usp(:,2))/0.511;

betae=zeros(nfiles,2,'double');
betae(:,1)=sqrt(1-1./gamae(:,1).^2);
betae(:,2)=sqrt(1-1./gamae(:,2).^2);

% Calculation of cooling force effect on ion beam
betaion0=zeros(nfiles,2,'double');
betaion0(:,1)=sqrt(1-1/gamaion^2);
betaion0(:,2)=sqrt(1-1/gamaion^2);
gamaion0=zeros(nfiles,2,'double');
gamaion0(:,1)=gamaion;
gamaion0(:,2)=gamaion;

% The relative energy
Erel=zeros(nfiles,2,'double');
Erel(:,1)=(sqrt(me^2*c^4+mion^2*c^4+2*me*mion*c^4*gamaion0(:,1).*gamae(:,1).*(1-betaion0(:,1).*betae(:,1)))-me*c^2-mion*c^2)/e;
Erel(:,2)=(sqrt(me^2*c^4+mion^2*c^4+2*me*mion*c^4*gamaion0(:,2).*gamae(:,2).*(1-betaion0(:,2).*betae(:,2)))-me*c^2-mion*c^2)/e;

%===========================================================


clear a b num 
%% Predefined variables
Ap1=zeros(nfiles,1,'double');
Ac1=zeros(nfiles,1,'double');
total=zeros(nfiles,2,'double');

Ap2=zeros(nfiles,1,'double');
Ac2=zeros(nfiles,1,'double');

Cell_DR_Counts_DCCT_Data=cell(nfiles,4);
% Read DR Data from local files (Each DR file consists of 8 columns)
% the 1st column stands for beam current of main ion beam (DCCT), the 2nd 
% column stands for Recombination ion counts (Counts), 
for i=1:nfiles   
    [DCCT,Count,NIM3,NIM2,NIM1,NIM5,NIM4,NIM6,Bump,MCPCount]=textread(filename{i},...
           '%f%f%f%f%f%f%f%f%f%f','headerlines',0); 

% Get the Size of DCCT
[m0,n0]=size(DCCT);

% Build zero matrix
index0=zeros(m0,1,'double');
index00=zeros(m0,1,'double');
index1=zeros(m0,1,'double');
index11=zeros(m0,1,'double');

    
% Judge where we start picking up data
for j=1:m0
    
    if j+cyclenum1<=m0 && NIM1(j,1)==1 
        for p=1:nummax
          if NIM2(p+j,1)~=0              
             index0(j,1)=1;
             index00(p+j,1)=1;             
          end
        end
    end
    
    if j+cyclenum1<=m0 && NIM3(j,1)==1 
        for k=1:nummax
          if NIM4(k+j,1)~=0
             index1(j,1)=1;
             index11(k+j,1)=1;             
          end
        end
    end
    
end

clear p k j 
% Get the points where DR data picked up
dex0=find(index0(:,1)==1);
dex00=find(index00(:,1)==1);

dex1=find(index1(:,1)==1);
dex11=find(index11(:,1)==1);

[a0,b0]=size(dex0);
[a1,b1]=size(dex1);

clear b0 b1
% if a0==a1
% Build DCCT and Counts matrix, DCCTs and Counts, respectively.
DCCTs1=zeros(a0*deltan,1,'double');
Counts1=zeros(a0*deltan,1,'double');
CoolCounts1=zeros(a0*deltamn1,1,'double');
CoolDCCTs1=zeros(a0*deltamn1,1,'double');

DCCTs2=zeros(a1*deltan,1,'double');
Counts2=zeros(a1*deltan,1,'double');
CoolCounts2=zeros(a1*deltamn1,1,'double');
CoolDCCTs2=zeros(a1*deltamn1,1,'double');    
    
% Select the real events of recombined ions by up voltages
for k=1:a0-1
    % Select the real counts of recombination
    Counts1((k-1)*deltan+1:(k-1)*deltan+deltan,1)=Count(dex0(k,1)+1+pick:dex0(k,1)+deltan+pick,1);
    % Select the real cooling counts of recombination
    CoolCounts1((k-1)*deltamn1+1:(k-1)*deltamn1+deltamn1,1)=Count(dex0(k+1,1)-deltamn1:dex0(k+1,1)-1,1);

    % Select the real DCCT of main beam ions       
if DCCT(dex0(k,1)+1,1)~=0
   DCCTs1((k-1)*deltan+1:(k-1)*deltan+deltan,1)=DCCT(dex0(k,1)+1,1);
           
elseif DCCT(dex0(k,1)+1,1)==0
    for m=1:cycletime1
        if (dex0(k,1)+1-m>0) && (dex0(k,1)+1+m<=m0) && (DCCT(dex0(k,1)+1+m,1) > DCCT(dex0(k,1)+1-m,1))
           DCCTs1((k-1)*deltan+1:(k-1)*deltan+deltan,1)=DCCT(dex0(k,1)+1+m,1);
        elseif (dex0(k,1)+1-m>0) && (dex0(k,1)+1+m<=m0) && (DCCT(dex0(k,1)+1+m,1) < DCCT(dex0(k,1)+1-m,1))
           DCCTs1((k-1)*deltan+1:(k-1)*deltan+deltan,1)=DCCT(dex0(k,1)+1-m,1);              
        end
    end
end

     
% Select the real Cool DCCT of primary beam ions      
if DCCT(dex00(k,1)+1,1)~=0
   CoolDCCTs1((k-1)*deltamn1+1:(k-1)*deltamn1+deltamn1,1)=DCCT(dex00(k,1)+1,1);
           
elseif DCCT(dex00(k,1)+1,1)==0
       for h=1:cycletime1
          if (dex00(k,1)+1-h>0) && (dex00(k,1)+1+h<=m0) && (DCCT(dex00(k,1)+1+h,1) > DCCT(dex00(k,1)+1-h,1))
             CoolDCCTs1((k-1)*deltamn1+1:(k-1)*deltamn1+deltamn1,1)=DCCT(dex00(k,1)+1+h,1);   
          elseif (dex00(k,1)+1-h>0) && (dex00(k,1)+1+h<=m0) && (DCCT(dex00(k,1)+1+h,1) < DCCT(dex00(k,1)+1-h,1))
             CoolDCCTs1((k-1)*deltamn1+1:(k-1)*deltamn1+deltamn1,1)=DCCT(dex00(k,1)+1-h,1);   
          end
       end
end

end


clear k h m
% Select the real events of recombined ions by down voltages
for k=1:a1-1
    % Select the real counts of recombination
    Counts2((k-1)*deltan+1:(k-1)*deltan+deltan,1)=Count(dex1(k,1)+1+pick:dex1(k,1)+deltan+pick,1);
    % Select the real cooling counts of recombination
    CoolCounts2((k-1)*deltamn1+1:(k-1)*deltamn1+deltamn1,1)=Count(dex1(k+1,1)-deltamn1:dex1(k+1,1)-1,1);

    % Select the real DCCT of main beam ions       
if DCCT(dex1(k,1)+1,1)~=0
   DCCTs2((k-1)*deltan+1:(k-1)*deltan+deltan,1)=DCCT(dex1(k,1)+1,1);
           
elseif DCCT(dex1(k,1)+1,1)==0
    for m=1:cycletime1
        if (dex1(k,1)+1-m>0) && (dex1(k,1)+1+m<=m0) && (DCCT(dex1(k,1)+1+m,1) > DCCT(dex1(k,1)+1-m,1))
           DCCTs2((k-1)*deltan+1:(k-1)*deltan+deltan,1)=DCCT(dex1(k,1)+1+m,1);
        elseif (dex1(k,1)+1-m>0) && (dex1(k,1)+1+m<=m0) && (DCCT(dex1(k,1)+1+m,1) < DCCT(dex1(k,1)+1-m,1))
           DCCTs2((k-1)*deltan+1:(k-1)*deltan+deltan,1)=DCCT(dex1(k,1)+1-m,1);              
        end
    end
end

     
% Select the real Cool DCCT of main beam ions      
if DCCT(dex11(k,1)+1,1)~=0
   CoolDCCTs2((k-1)*deltamn1+1:(k-1)*deltamn1+deltamn1,1)=DCCT(dex11(k,1)+1,1);
           
elseif DCCT(dex11(k,1)+1,1)==0
       for h=1:cycletime1
          if (dex11(k,1)+1-h>0) && (dex11(k,1)+1+h<=m0) && (DCCT(dex11(k,1)+1+h,1) > DCCT(dex11(k,1)+1-h,1))
             CoolDCCTs2((k-1)*deltamn1+1:(k-1)*deltamn1+deltamn1,1)=DCCT(dex11(k,1)+1+h,1);   
          elseif (dex11(k,1)+1-h>0) && (dex11(k,1)+1+h<=m0) && (DCCT(dex11(k,1)+1+h,1) < DCCT(dex11(k,1)+1-h,1))
             CoolDCCTs2((k-1)*deltamn1+1:(k-1)*deltamn1+deltamn1,1)=DCCT(dex11(k,1)+1-h,1);   
          end
       end
end

end

clear k h m
%else
%       display('Search of detuning data mistakes at the data file of');
%       display(filename{i});

% end

% Find the bad data of Labview or data acquisition card created in the
% sample process
% Select the real detuning data meet the wanted conditions by up voltages
Bad1Counts1=find(Counts1(:,1)>900);
if isempty(Bad1Counts1)==0
Counts1(Bad1Counts1)=[];
DCCTs1(Bad1Counts1)=[];
end

Bad1DCCTs1=find(DCCTs1(:,1)<5);
if isempty(Bad1DCCTs1)==0
Counts1(Bad1DCCTs1)=[];
DCCTs1(Bad1DCCTs1)=[];  
end


% Select the real cool data meet the wanted conditions
Bad1CoolCounts1=find(CoolCounts1(:,1)>900);
if isempty(Bad1CoolCounts1)==0
CoolCounts1(Bad1CoolCounts1)=[];
CoolDCCTs1(Bad1CoolCounts1)=[];
end

Bad2CoolDCCTs1=find(CoolDCCTs1(:,1)<5);
if isempty(Bad2CoolDCCTs1)==0
CoolCounts1(Bad2CoolDCCTs1)=[];
CoolDCCTs1(Bad2CoolDCCTs1)=[];
end


% Select the real detuning data meet the wanted conditions by down voltages
Bad1Counts2=find(Counts2(:,1)>900);
if isempty(Bad1Counts2)==0
Counts2(Bad1Counts2)=[];
DCCTs2(Bad1Counts2)=[];
end

Bad2DCCTs2=find(DCCTs2(:,1)<5);
if isempty(Bad2DCCTs2)==0
Counts2(Bad2DCCTs2)=[];
DCCTs2(Bad2DCCTs2)=[];  
end

% Select the real cool data meet the wanted conditions
Bad1CoolCounts2=find(CoolCounts2(:,1)>900);
if isempty(Bad1CoolCounts2)==0
CoolCounts2(Bad1CoolCounts2)=[];
CoolDCCTs2(Bad1CoolCounts2)=[];
end

Bad2CoolDCCTs2=find(CoolDCCTs2(:,1)<5);
if isempty(Bad2CoolDCCTs2)==0
CoolCounts2(Bad2CoolDCCTs2)=[];
CoolDCCTs2(Bad2CoolDCCTs2)=[];
end

% Save all the Counts data in a cell array {Cell_DR_Counts_DCCT_Data}
Cell_DR_Counts_DCCT_Data{i,1}(:,2)=Counts1(:,1);
Cell_DR_Counts_DCCT_Data{i,1}(:,1)=DCCTs1(:,1);
Cell_DR_Counts_DCCT_Data{i,2}(:,2)=Counts2(:,1);
Cell_DR_Counts_DCCT_Data{i,2}(:,1)=DCCTs2(:,1);

Cell_DR_Counts_DCCT_Data{i,3}(:,2)=CoolCounts1(:,1);
Cell_DR_Counts_DCCT_Data{i,3}(:,1)=CoolDCCTs1(:,1);
Cell_DR_Counts_DCCT_Data{i,4}(:,2)=CoolCounts2(:,1);
Cell_DR_Counts_DCCT_Data{i,4}(:,1)=CoolDCCTs2(:,1);

% normalnize the counts by DCCT
Ap1(i,1)=mean(Counts1(:,1)./DCCTs1(:,1));
[m1,n1]=size(Counts1);
total(i,1)=mean(m1*Counts1(:,1));

Ac1(i,1)=mean(CoolCounts1(:,1)./CoolDCCTs1(:,1));

Ap2(i,1)=mean(Counts2(:,1)./DCCTs2(:,1));
[m2,n2]=size(Counts2);
total(i,2)=mean(m2*Counts2(:,1));

Ac2(i,1)=mean(CoolCounts2(:,1)./CoolDCCTs2(:,1));

clear m1 n1 m2 n2 m11 n11 m22 n22

end


%% Relative rate coefficients
Alpha_Up=zeros(nfiles,2,'double');
Alpha_Down=zeros(nfiles,2,'double');

Alpha_Up(:,1)=Ap1(:,1).*betae(:,1).*betaion0(:,1).*gamaion0(:,1).^2;
Alpha_Up(:,2)=Ac1(:,1)*betaion*betaion*gamaion^2;

Alpha_Down(:,1)=Ap2(:,1).*betae(:,2).*betaion0(:,2).*gamaion0(:,2).^2;
Alpha_Down(:,2)=Ac2(:,1)*betaion*betaion*gamaion^2;

Alpha0=zeros(nfiles,2,'double');
if type==1
    Alpha0(:,1)=Alpha_Up(:,1).*qcharge*e^2*c^2*pi*(radius_e/100)^2*10^15/Ie/lc;      %in unit of cm^3/s 
    Alpha0(:,2)=Alpha_Down(:,1).*qcharge*e^2*c^2*pi*(radius_e/100)^2*10^15/Ie/lc;
else    
    Alpha0(:,1)=Alpha_Up(:,1)./Alpha_Up(:,2);
    Alpha0(:,2)=Alpha_Down(:,1)./Alpha_Down(:,2);
end

Alpha0_All=[detuning_voltages,total,Erel,Alpha0];
Alpha_All=[Alpha_Up,Alpha_Down];
%% Plots
figure(1)
plot(Erel(:,1),Alpha0(:,1));
title('Alpha0 of up detuning')

figure(2)
plot(Erel(:,2),Alpha0(:,2));
title('Alpha0 of down detuning')

% toc; % Obtain the elapsed time of this program

end