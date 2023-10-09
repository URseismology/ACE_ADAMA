%%%%%%% Fig suppliment: plot the boxplots of map velocity (by depth) for sequencer/K-Mean class results
%%%%%%% Siyu Xue -- SEP 5. 2023
clear

%%%% INPUT the depth of interest
mdepth = 20;    % [km] in depth

%%% NOTE: the results are not so ideal, so probably won't use this code

%% (1.1) Prepare the ADAMA Velocity map (slice)

load('/Users/sxue3/Documents/BayMap_Figures/Data/Litho1D.mat')  % load Litho 1D models
load('/Users/sxue3/Documents/BayMap_Figures/Data/ADAMAinvert/ADAMA_1D_fixed1p7.mat')  % load ADAMA inversion results
load('/Users/sxue3/Documents/BayMap_Figures/Data/AfrDownSampleDots.mat')  % load data locations
load('/Users/sxue3/Documents/BayMap_Figures/Data/RayDataQly.mat')   % load ADAMA data quality
% % load the coast
% afcoast = load('./Data/AfricaCoast.mat');
% aflon = wrapTo180(afcoast.XY(:, 1));
% aflat = afcoast.XY(:,2);
% coastline = polyshape(aflon, aflat);
% % load African country boundaries
% S1 = load('./Data/AfrCountry.mat').S1;
% madgline = polyshape((S1(448).X).',(S1(448).Y).');  % load Madgascar's Big Island
% 
% % load cratons
% load('/Users/sxue3/Documents/BayMap_Figures/Data/AfricanCratons.mat');
% % Sahara
% Sahara_line = polyshape(AC.Sahara(:,2), AC.Sahara(:,1));
% 
% % Congo
% Congo_line = polyshape(AC.Congo(:,2), AC.Congo(:,1));
% 
% % Tanza
% Tanza_line = polyshape(AC.Tanza(:,2), AC.Tanza(:,1));
% 
% % Kala
% Kala_line = polyshape(AC.Kala(:,2), AC.Kala(:,1));
% 
% % West
% West_line = polyshape(AC.West(:,2), AC.West(:,1));

%% (1.2) Read in the classification results 
% note : the classofication has the same ordering as dslat, dslon, Litho & ADAMA vels
C1model = readmatrix('/Users/sxue3/Documents/BayMap_Figures/Data/Sequencer/AfricaClass_k4.csv');
C1model = C1model + 1;

%% (2.1) Assemble Litho & ADAMA slices

iD = mdepth*2+1;    % depth = 500*(iD-1) [m]


% Assemble Litho slice
vel2D_litho = zeros(size(dslon));

for ip = 1:size(dslon,1) 
    Lithovs = Litho1D{ip}.pvs;
    vel2D_litho(ip) = Lithovs(iD);
end


% Assemble ADAMA slice
vel2D_ADAMA = zeros(size(dslon));

for ip = 1:size(dslon,1) 
    if ~isempty(ADAMA1D{ip})  % use ADAMA if there is an update
        ADAMAvs = ADAMA1D{ip}.vsv_update;
        fupdate = find(ADAMAvs(:,1) ~= 0);
        vel2D_ADAMA(ip) = ADAMAvs(fupdate(end),iD);
    else
        vel2D_ADAMA(ip) = vel2D_litho(ip);
    end
end


%% (3.1) Plot the crustal velocity boxplot using Crust1.0

figure(310)
subplot(1,2,1)
% plot the velocity distribution of ADAMA classes
histogram(vel2D_ADAMA(C1model==1)/1000, 'FaceAlpha',0.5, 'NumBins',100, 'EdgeAlpha', 0, 'FaceColor',[0 0.447 0.741])
hold on
histogram(vel2D_ADAMA(C1model==2)/1000, 'FaceAlpha',0.5, 'NumBins',100, 'EdgeAlpha', 0, 'FaceColor',[0.581 0.325 0.098])
histogram(vel2D_ADAMA(C1model==3)/1000, 'FaceAlpha',0.5, 'NumBins',100, 'EdgeAlpha', 0, 'FaceColor',[0.929 0.694 0.125])
histogram(vel2D_ADAMA(C1model==4)/1000, 'FaceAlpha',0.5, 'NumBins',100, 'EdgeAlpha', 0, 'FaceColor',[0.8 0.8 0.8])
title(['ADAMA Dist: k=4 at ', num2str(mdepth), 'km'])
legend('C1', 'C2','C3', 'C4')
xlabel('Velocity (km/s)')
xlim([3.5 4.5])

Litho_updates = (vel2D_ADAMA-vel2D_litho)./vel2D_litho*100;
subplot(1,2,2)
% plot the distribution of 
histogram(Litho_updates(C1model==1)/1000, 'FaceAlpha',0.5, 'NumBins',100, 'EdgeAlpha', 0)
hold on
histogram(Litho_updates(C1model==2)/1000, 'FaceAlpha',0.5, 'NumBins',100, 'EdgeAlpha', 0)
histogram(Litho_updates(C1model==3)/1000, 'FaceAlpha',0.5, 'NumBins',100, 'EdgeAlpha', 0)
histogram(Litho_updates(C1model==4)/1000, 'FaceAlpha',0.5, 'NumBins',100, 'EdgeAlpha', 0)
title(['Litho Update %: k=4 at ', num2str(mdepth), 'km'])
xlim([-0.05, 0.05])

% set the size of the plot
x0=300;
y0=200;
fwidth=800;
fheight=350;
set(gcf,'position',[x0,y0,fwidth,fheight]);

%% statistics

class1_mean = mean(vel2D_ADAMA(C1model == 1))/1000;
class2_mean = mean(vel2D_ADAMA(C1model == 2))/1000;
class3_mean = mean(vel2D_ADAMA(C1model == 3))/1000;
class4_mean = mean(vel2D_ADAMA(C1model == 4))/1000;

class1_med = median(vel2D_ADAMA(C1model == 1))/1000;
class2_med = median(vel2D_ADAMA(C1model == 2))/1000;
class3_med = median(vel2D_ADAMA(C1model == 3))/1000;
class4_med = median(vel2D_ADAMA(C1model == 4))/1000;

class1_var = var(vel2D_ADAMA(C1model == 1)/1000);
class2_var = var(vel2D_ADAMA(C1model == 2)/1000);
class3_var = var(vel2D_ADAMA(C1model == 3)/1000);
class4_var = var(vel2D_ADAMA(C1model == 4)/1000);






