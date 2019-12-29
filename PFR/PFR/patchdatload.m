function [Data,Adds]=patchdatload(Start,Incr,End,Datnum)
% Code to deal with opening Data files that follow a (data_Iteration#_Data
% File) scheme
Bark = (Start:Incr:End)';
[r,~] = size(Bark);
% names=zeros(r);
names=cell(size(1:r));

j=1; % Counter to move through the name array
for i = Start:Incr:End
    
    string=strcat('results/data_',num2str(i),'_',num2str(Datnum),'.txt');
    names{j}=string;
    j=j+1;
end

% Reads the first name array and grabs the experimental values
Temp=dlmread(names{1});
Data(:,1)=Temp(:,4);

% Reads the rest of the arrays and grabs the simulated values
for i=1:(j-1)
    Temp=dlmread(names{i});
    Data(:,i+1) = Temp(:,5);
end

% Loading Data into MATLAB from Text Files
% Ft='.txt'; % Designates File Type
% FileTemp=num2str(FileTemp); % Turns the Temp Value into a String
% DATA='data'; % Data string
% InDat='data/Hauchhum_et_Mahanta_13XZeolite_';
% FullDatString=[PATH DATA Time '_' FileTemp Ft];
% DatString='data/Hauchhum_et_Mahanta_13XZeolite_';
% DatString=[DatString num2str(FileTemp) 'C' Ft];
% DIN=dlmread(DatString,'',1,0);

% Info from Data File
% T=DIN(1,3); % Grabs the temp [C]
% T=Temp(1,3); % Grabs the temp [C]
% T=T+273.15; % Temp in [K]
% Pco2=DIN(:,1); % Grabs CO2 Pressure in [Pa]
Pco2=Temp(:,1); % Grabs CO2 Pressure [Pa]
% EWF=DIN(:,4); % Grabs Experimental Weight Fraction
% EWF=Temp(:,4); % Grabs Experimental Weight Fraction
% Adds.Temp=T;
Adds.PCO2=Pco2;
% Adds.ExpWeightFrac=EWF;
end