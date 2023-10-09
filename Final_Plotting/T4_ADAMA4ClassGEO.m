%%%%%% Compute statistics of ADAMA velocity for each class, given GEO classification 
%%%%%% 
clear

%% (1) Read in data
load('/Users/sxue3/Documents/BayMap_Figures/Data/Litho1D.mat')  % load Litho 1D models
load('/Users/sxue3/Documents/BayMap_Figures/Data/ADAMAinvert/ADAMA_1D_fixed1p7.mat')  % load ADAMA inversion results
load('/Users/sxue3/Documents/BayMap_Figures/Data/AfrDownSampleDots.mat')  % load data locations

% load the African GEOs
load('/Users/sxue3/Documents/BayMap_Figures/Data/GeoData/AfricanCrustal.mat')

% load the ADAMA classes from KMean
Aclass = readmatrix('/Users/sxue3/Documents/BayMap_Figures/Data/Sequencer/AfricaClass_k4_ds1.csv');


% load the African coast boundaries
afcoast = load('./Data/GeoData/AfricaCoast.mat');
aflon = wrapTo180(afcoast.XY(:, 1));
aflat = afcoast.XY(:,2);
coastline = polyshape(aflon, aflat);
S1 = load('./Data/GeoData/AfrCountry.mat').S1;
madgline = polyshape((S1(448).X).',(S1(448).Y).');  % load Madgascar's Big Island

% load the ADAMA Vel classes
Aclass = readmatrix('/Users/sxue3/Documents/BayMap_Figures/Data/Sequencer/AfricaClass_k4.csv');

%% (2) Get ADAMA Data

vel1D_ADAMA = zeros(length(dslon), 4); % store all vel to depth=50km with 0.5km resolution

for ip = 1:size(dslon,1) 
    if ~isempty(ADAMA1D{ip})  % use ADAMA if there is an update
        ADAMAvs = ADAMA1D{ip}.vsv_update;
        fupdate = find(ADAMAvs(:,1) ~= 0);
        vel1D_ADAMA(ip, :) = ADAMAvs(fupdate(end),[21, 41, 61, 81])/1000;
    else
        Lithovs = Litho1D{ip}.pvs;
        vel1D_ADAMA(ip, :) = Lithovs([21, 41, 61, 81])/1000;
    end
end


%% (3) Archean Blocks Stats

LeoManShield = polyshape(GEO.Archean.WestAfrica.LeoManShield(:,2), GEO.Archean.WestAfrica.LeoManShield(:,1));
Reguibatshield = polyshape(GEO.Archean.WestAfrica.Reguibatshield(:,2), GEO.Archean.WestAfrica.Reguibatshield(:,1));

in_Archean_WA = (isinterior(LeoManShield, [dslon, dslat]) | isinterior(Reguibatshield, [dslon, dslat]));

disp(['% in Africa: ', num2str(sum(in_Archean_WA)/2381*100)]);      % % of Archean_WA in Africa

temp_class = in_Archean_WA .* Aclass;

area_C1 = sum(temp_class == 1)/sum(in_Archean_WA)*100;
area_C2 = sum(temp_class == 2)/sum(in_Archean_WA)*100;
area_C3 = sum(temp_class == 3)/sum(in_Archean_WA)*100;
area_C4 = sum(temp_class == 4)/sum(in_Archean_WA)*100;
disp(['% of Class 1: ', num2str(area_C1)]);
disp(['% of Class 2: ', num2str(area_C2)]);
disp(['% of Class 3: ', num2str(area_C3)]);
disp(['% of Class 4: ', num2str(area_C4)]);

C1_avg = mean(vel1D_ADAMA(temp_class == 1, :));
C2_avg = mean(vel1D_ADAMA(temp_class == 2, :));
C3_avg = mean(vel1D_ADAMA(temp_class == 3, :));
C4_avg = mean(vel1D_ADAMA(temp_class == 4, :));
disp('C1 average for depth 10, 20, 30, 40km: ')
disp(['     ', num2str(C1_avg)])
disp('C2 average for depth 10, 20, 30, 40km: ')
disp(['     ', num2str(C2_avg)])
disp('C3 average for depth 10, 20, 30, 40km: ')
disp(['     ', num2str(C3_avg)])
disp('C3 average for depth 10, 20, 30, 40km: ')
disp(['     ', num2str(C4_avg)])

% special case since there is only 1 DP for Class3
D1_avg = (C1_avg(1)*area_C1+C2_avg(1)*area_C2+3.6211*area_C3+C4_avg(1)*area_C4)/100;
D2_avg = (C1_avg(2)*area_C1+C2_avg(2)*area_C2+3.8287*area_C3+C4_avg(2)*area_C4)/100;
D3_avg = (C1_avg(3)*area_C1+C2_avg(3)*area_C2+4.0981*area_C3+C4_avg(3)*area_C4)/100;
D4_avg = (C1_avg(4)*area_C1+C2_avg(4)*area_C2+4.7432*area_C3+C4_avg(4)*area_C4)/100;
disp(['Weighted average for Class 1: ', num2str(D1_avg)])
disp(['Weighted average for Class 2: ', num2str(D2_avg)])
disp(['Weighted average for Class 3: ', num2str(D3_avg)])
disp(['Weighted average for Class 4: ', num2str(D4_avg)])





%%
% %%%% For Checking.....
% hold on
% % plot continent
% plot(coastline, 'EdgeColor', 'k', 'LineWidth', 4, 'FaceAlpha', 0);  % plot the African boundary
% plot(madgline, 'EdgeColor', 'k', 'LineWidth', 4, 'FaceAlpha', 0);  % plot the African boundary
% 
% hold on
% % plot continent
% plot(LeoManShield, 'EdgeColor', 'r', 'LineWidth', 4, 'FaceAlpha', 0);  % plot the African boundary
% plot(Reguibatshield, 'EdgeColor', 'b', 'LineWidth', 4, 'FaceAlpha', 0);  % plot the African boundary
% 
% plot(dslon(in_Archean_WA), dslat(in_Archean_WA), 'g*')

%% (2) get the locations inside each region


% %% save geo data from CVS to mat
% 
% Cap = readmatrix('/Users/sxue3/Documents/BayMap_Figures/Data/GeoData/Mobile_Belts/capbelt.csv');
% Kibaran = readmatrix('/Users/sxue3/Documents/BayMap_Figures/Data/GeoData/Mobile_Belts/Kibaranbelt.csv');
% Namaquanatal = readmatrix('/Users/sxue3/Documents/BayMap_Figures/Data/GeoData/Mobile_Belts/Namaquanatalbelt.csv');
% Ruwenzky = readmatrix('/Users/sxue3/Documents/BayMap_Figures/Data/GeoData/Mobile_Belts/ruwenzkybelt.csv');
% Mauritania = readmatrix('/Users/sxue3/Documents/BayMap_Figures/Data/GeoData/Mobile_Belts/Mauritaniabelt.csv');
% 
% GEO.MobileBelts.Other.Cap_belt = Cap;
% GEO.MobileBelts.Other.Kibaran_belt = Kibaran;
% GEO.MobileBelts.Other.Namaquanatal = Namaquanatal;
% GEO.MobileBelts.Other.Ruwenzky = Ruwenzky;
% GEO.MobileBelts.Other.Mauritania = Mauritania;
% 
% GEO.MobileBelts.WestAfrica_belt = readmatrix('/Users/sxue3/Documents/BayMap_Figures/Data/GeoData/Mobile_Belts/WestAfricaBelt.csv');
% GEO.MobileBelts.Ouban_Damara.Damara = readmatrix('/Users/sxue3/Documents/BayMap_Figures/Data/GeoData/Mobile_Belts/Damarabelt.csv');
% GEO.MobileBelts.Ouban_Damara.Oubangides = readmatrix('/Users/sxue3/Documents/BayMap_Figures/Data/GeoData/Mobile_Belts/Oubangidesbelt.csv');
% 
% GEO.Columns = 'The first column is latitude, and the second column is longtitude';
% 
% 
% GEO.Orogens.AtlasMountains = readmatrix('/Users/sxue3/Documents/BayMap_Figures/Data/GeoData/Orogens/Atlas.csv');
% GEO.Orogens.EastAfrica.EastAfrica = readmatrix('/Users/sxue3/Documents/BayMap_Figures/Data/GeoData/Orogens/East_Afr_Orogenic_Zone.csv');
% GEO.Orogens.EastAfrica.NorthMozamb = readmatrix('/Users/sxue3/Documents/BayMap_Figures/Data/GeoData/Orogens/North_Mozamb_Belt.csv');
% 
% 
% GEO.Basins.Congo = readmatrix('/Users/sxue3/Documents/BayMap_Figures/Data/GeoData/Basins/congo_basin.csv');
% 
% GEO.Basins.Tindouf.Tindouf = readmatrix('/Users/sxue3/Documents/BayMap_Figures/Data/GeoData/Basins/Tindouf_basin.csv');
% GEO.Basins.Tindouf.MauriSenba = readmatrix('/Users/sxue3/Documents/BayMap_Figures/Data/GeoData/Basins/MauritaniaSenbasin.csv');
% 
% GEO.Basins.Taoudeni.B1 = readmatrix('/Users/sxue3/Documents/BayMap_Figures/Data/GeoData/Basins/TaoudeniBasin.csv');
% GEO.Basins.Taoudeni.B2 = readmatrix('/Users/sxue3/Documents/BayMap_Figures/Data/GeoData/Basins/Taoudeni_basin2.csv');
% 
% 
% GEO.Archean.WestAfrica.LeoManShield = readmatrix('/Users/sxue3/Documents/BayMap_Figures/Data/GeoData/Archean_Blocks/LeoManShield_WA1.csv');
% GEO.Archean.WestAfrica.Reguibatshield = readmatrix('/Users/sxue3/Documents/BayMap_Figures/Data/GeoData/Archean_Blocks/Reguibatshield_WA2.csv');
% 
% GEO.Archean.Tanzania = readmatrix('/Users/sxue3/Documents/BayMap_Figures/Data/GeoData/Archean_Blocks/Tanzaniacraton.csv');
% 
% GEO.Archean.CongoAll.AC = readmatrix('/Users/sxue3/Documents/BayMap_Figures/Data/GeoData/Archean_Blocks/AC.csv');
% GEO.Archean.CongoAll.Bomu = readmatrix('/Users/sxue3/Documents/BayMap_Figures/Data/GeoData/Archean_Blocks/Bomu_craton.csv');
% GEO.Archean.CongoAll.Congo = readmatrix('/Users/sxue3/Documents/BayMap_Figures/Data/GeoData/Archean_Blocks/Congo_craton_3.csv');
% GEO.Archean.CongoAll.Kasai = readmatrix('/Users/sxue3/Documents/BayMap_Figures/Data/GeoData/Archean_Blocks/Kasai_craton.csv');
% GEO.Archean.CongoAll.Ntem = readmatrix('/Users/sxue3/Documents/BayMap_Figures/Data/GeoData/Archean_Blocks/Ntem_craton.csv');
% 
% GEO.Archean.SaharaMeta.A1 = readmatrix('/Users/sxue3/Documents/BayMap_Figures/Data/GeoData/Archean_Blocks/SMC_Archeanblock1.csv');
% GEO.Archean.SaharaMeta.A2 = readmatrix('/Users/sxue3/Documents/BayMap_Figures/Data/GeoData/Archean_Blocks/SMCArcheanblock2.csv');
% GEO.Archean.SaharaMeta.A3 = readmatrix('/Users/sxue3/Documents/BayMap_Figures/Data/GeoData/Archean_Blocks/SMCArcheanblock3.csv');
% GEO.Archean.SaharaMeta.A4 = readmatrix('/Users/sxue3/Documents/BayMap_Figures/Data/GeoData/Archean_Blocks/SMCArcheanblock4.csv');
% 
% GEO.Archean.WestAfrica_MZ.Dahomeyshield = readmatrix('/Users/sxue3/Documents/BayMap_Figures/Data/GeoData/Archean_Blocks/Dahomeyshield_MZ.csv');
% GEO.Archean.WestAfrica_MZ.Tuaregshield = readmatrix('/Users/sxue3/Documents/BayMap_Figures/Data/GeoData/Archean_Blocks/Tuaregshield_MZ.csv');
% GEO.Archean.KaapVaal = readmatrix('/Users/sxue3/Documents/BayMap_Figures/Data/GeoData/Archean_Blocks/Zimb_Kaapvaal_Craton.csv');
% 
