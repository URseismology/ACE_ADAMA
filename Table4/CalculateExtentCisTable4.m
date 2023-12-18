%%%
%  Table 4 --- Area Norm -- Extent %
%   C1,C2,C3,C4
%   AB: Archean Blocks
%   B: Basins
%   MB: Mobile Belts
%   O: Orogens
%   U: Unclassified

% J-J legre
%%

clear 
clc 

Data = table2array(readtable('DataAfricaTable4.csv'));
Data(isnan(Data))=0; % replace nan values

% Arrays for each Geologic group
DataAB = Data(2:7,1:25);
DataB = Data(9:11,1:25);
DataMB = Data(13:15,1:25);
DataO = Data(17:18,1:25);
DataU = Data(20:22,1:25);

GroupCs_Extent = zeros(6,4);
DataGroups = vertcat(DataAB,DataB,DataMB,DataO,DataU);


%%
for e=2:5  
    
    %i[1,4]
    Sum_Ci_AB = sum((DataAB(:,1) .*  DataAB(:,e))) ;
    Sum_Ci_All = sum((DataGroups(:,1) .*  DataGroups(:,e)));
    
    %AB
    Ci_AB = round((Sum_Ci_AB / Sum_Ci_All)*100,1);
    GroupCs_Extent(1,e-1) = Ci_AB ;
    
    %B
    Sum_Ci_B = sum((DataB(:,1) .*  DataB(:,e))) ;
    Ci_B = round((Sum_Ci_B / Sum_Ci_All)*100,1);
    GroupCs_Extent(2,e-1) = Ci_B; 
    
    %MB
    Sum_Ci_MB = sum((DataMB(:,1) .*  DataMB(:,e))) ;
    Ci_MB = round((Sum_Ci_MB / Sum_Ci_All)*100,1);
    GroupCs_Extent(3,e-1) = Ci_MB ;
    
    %O
    Sum_Ci_O = sum((DataO(:,1) .*  DataO(:,e))) ;
    Ci_O = round((Sum_Ci_O / Sum_Ci_All)*100,1);
    GroupCs_Extent(4,e-1) = Ci_O ;
    
    %U
    Sum_Ci_U = sum((DataU(:,1) .*  DataU(:,e))) ;
    Ci_U = round((Sum_Ci_U / Sum_Ci_All)*100,1);
    GroupCs_Extent(5,e-1) = Ci_U ;
    
    %Ci
    Total_Ci = Ci_AB + Ci_B + Ci_MB + Ci_O + Ci_U; %sum should be 100
    GroupCs_Extent(6,e-1) = Total_Ci; %Vertical
             
end 

GroupCs_Extent


