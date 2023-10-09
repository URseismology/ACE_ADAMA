%%%%%% Prepare the data for Sequencer and plot the Sequencer index results
clear

%% (1) Read in data
load('/Users/sxue3/Documents/BayMap_Figures/Data/Litho1D.mat')  % load Litho 1D models
load('/Users/sxue3/Documents/BayMap_Figures/Data/ADAMAinvert/ADAMA_1D_fixed1p7.mat')  % load ADAMA inversion results
load('/Users/sxue3/Documents/BayMap_Figures/Data/AfrDownSampleDots.mat')  % load data locations
load('/Users/sxue3/Documents/BayMap_Figures/Data/GeoData/AfricanCrustal.mat') % load crustal types

% load the African coast boundaries
afcoast = load('./Data/GeoData/AfricaCoast.mat');
aflon = wrapTo180(afcoast.XY(:, 1));
aflat = afcoast.XY(:,2);
coastline = polyshape(aflon, aflat);
S1 = load('./Data/GeoData/AfrCountry.mat').S1;
madgline = polyshape((S1(448).X).',(S1(448).Y).');  % load Madgascar's Big Island

% load cratons
load('/Users/sxue3/Documents/BayMap_Figures/Data/GeoData/AfricanCratons.mat');

%% (2) Prepare ADAMA vel by depth to a matrix

vel1D_ADAMA = zeros(length(dslon), 101); % store all vel to depth=50km with 0.5km resolution

for ip = 1:size(dslon,1) 
    if ~isempty(ADAMA1D{ip})  % use ADAMA if there is an update
        ADAMAvs = ADAMA1D{ip}.vsv_update;
        fupdate = find(ADAMAvs(:,1) ~= 0);
        vel1D_ADAMA(ip, :) = ADAMAvs(fupdate(end),1:101);
    else
        Lithovs = Litho1D{ip}.pvs;
        vel1D_ADAMA(ip, :) = Lithovs(1:101);
    end
end

%% (3) Smaller vel matrices for each individual craton

% Sahara
Sahara_line = polyshape(AC.Sahara(:,2), AC.Sahara(:,1));
in_sahara = isinterior(Sahara_line, [dslon, dslat]);
vel_sahara = vel1D_ADAMA(in_sahara, :);

% Congo
Congo_line = polyshape(AC.Congo(:,2), AC.Congo(:,1));
in_congo = isinterior(Congo_line, [dslon, dslat]);
vel_congo = vel1D_ADAMA(in_congo, :);

% Tanza
Tanza_line = polyshape(AC.Tanza(:,2), AC.Tanza(:,1));
in_tanza = isinterior(Tanza_line, [dslon, dslat]);
vel_tanza = vel1D_ADAMA(in_tanza, :);

% Kala
Kala_line = polyshape(AC.Kala(:,2), AC.Kala(:,1));
in_kala = isinterior(Kala_line, [dslon, dslat]);
vel_kala = vel1D_ADAMA(in_kala, :);

% West
West_line = polyshape(AC.West(:,2), AC.West(:,1));
in_west = isinterior(West_line, [dslon, dslat]);
vel_west = vel1D_ADAMA(in_west, :);

%% (4) Save the vel matrices to mat

% save('For_Sequencer.mat', 'vel_west', 'vel_kala', 'vel_tanza','vel_congo','vel_sahara','dslon','dslat','vel1D_ADAMA')


%% (5.1) Create the resulting Sequencer idx map

% read in the sequencer idx
IDX = readmatrix('/Users/sxue3/Documents/BayMap_Figures/Data/Sequencer/final_sequence_dataset_mat_vel1D_ADAMA_Seq1.txt');

% re-order the index
[a, b] = sort(IDX(:,2));
vel_idx = IDX(b);

xgrid = -20:0.1:55;    % lon
ygrid = -35:0.1:40;    % lat
[xq, yq] = meshgrid(xgrid, ygrid);
F = scatteredInterpolant(dslon,dslat,vel_idx, 'linear', 'none');
Vq_idx = F(xq,yq);


in_africa = ones(size(Vq_idx));

for ix = 1:size(Vq_idx, 1)
    for iy = 1:size(Vq_idx, 1)
        if isinterior(coastline, [xgrid(ix), ygrid(iy)]) ||  isinterior(madgline, [xgrid(ix), ygrid(iy)])
            continue
        else
            in_africa(iy, ix) = NaN; 
        end
    end
end

Vq_idx = Vq_idx .* in_africa;

%% (5.2) Plot the map showing index 
pcolor(xq,yq, Vq_idx)
shading flat
c1 = colorbar;
c1.Title.String = 'Index';
colormap('turbo')

hold on
% plot continent
plot(coastline, 'EdgeColor', 'k', 'LineWidth', 4, 'FaceAlpha', 0);  % plot the African boundary
plot(madgline, 'EdgeColor', 'k', 'LineWidth', 4, 'FaceAlpha', 0);  % plot the African boundary
% plot African cratons
plot(AC.Congo(:,2), AC.Congo(:,1), 'color', 'k', 'lineWidth', 2);
plot(AC.Kala(:,2), AC.Kala(:,1), 'color', 'k', 'lineWidth', 2);
plot(AC.Sahara(:,2), AC.Sahara(:,1), 'color', 'k', 'lineWidth', 2);
plot(AC.Tanza(:,2), AC.Tanza(:,1), 'color', 'k', 'lineWidth', 2);
plot(AC.West(:,2), AC.West(:,1), 'color', 'k', 'lineWidth', 2);
plot(AC.Zim(:,2), AC.Zim(:,1), 'color', 'k', 'lineWidth', 2);

xlabel('Longitude')
ylabel('Latitude')
minlon = -20;
maxlon = 55;
minlat = -36;
maxlat = 40;
xlim([minlon maxlon]);
ylim([minlat maxlat]);
axesm('mercator','maplatlimit',[minlat maxlat],'maplonlimit',[minlon maxlon]);
xticks(10*ceil(minlon/10):10:10*floor(maxlon/10));
yticks(10*ceil(minlat/10):10:10*floor(maxlat/10));
title('Sequencer Index')
set(gca,'fontsize', 14) 
set(gcf, 'color', 'w')


%% (6.1) Downsample 1D velocity 

%%%%%% ds1
vel1D_ADAMA_ds = [vel1D_ADAMA(:,1:10), vel1D_ADAMA(:,15), vel1D_ADAMA(:,20), vel1D_ADAMA(:,25)...
    , vel1D_ADAMA(:,35), vel1D_ADAMA(:,45), vel1D_ADAMA(:,55), vel1D_ADAMA(:,65), vel1D_ADAMA(:,75)...
    , vel1D_ADAMA(:,85)];

%%%%%% ds2
% vel1D_ADAMA_ds = [vel1D_ADAMA(:,1:15), vel1D_ADAMA(:,20), vel1D_ADAMA(:,25)...
%     , vel1D_ADAMA(:,35), vel1D_ADAMA(:,45), vel1D_ADAMA(:,55), vel1D_ADAMA(:,65), vel1D_ADAMA(:,75)...
%     , vel1D_ADAMA(:,85)];

%%%%%% ds3
% vel1D_ADAMA_ds = [vel1D_ADAMA(:,1:20), vel1D_ADAMA(:,25)...
%     , vel1D_ADAMA(:,35), vel1D_ADAMA(:,45), vel1D_ADAMA(:,55), vel1D_ADAMA(:,65), vel1D_ADAMA(:,75)...
%     , vel1D_ADAMA(:,85)];

pcolor(vel1D_ADAMA_ds)
shading flat

writematrix(round(vel1D_ADAMA_ds),'vel1D_ADAMA_ds1.cvs','Delimiter',',', 'FileType' ,'text')

% plot and check the downsampling
subplot(1,2,1)
pcolor(vel1D_ADAMA(IDX(:,2)+1,:))
shading flat

subplot(1,2,2)
pcolor(vel1D_ADAMA_ds(IDX(:,2)+1,:))
shading flat



%% (7.1) Read in the Africa Classification (k-mean clustering)

% read in the sequencer idx
Afrclass = readmatrix('/Users/sxue3/Documents/BayMap_Figures/Data/Sequencer/AfricaClass_k4.csv');

% re-order the index
xgrid = -20:0.1:55;    % lon
ygrid = -35:0.1:40;    % lat
[xq, yq] = meshgrid(xgrid, ygrid);
F = scatteredInterpolant(dslon,dslat,Afrclass, 'nearest', 'none');
Vq_class = F(xq,yq);

Vq_class = Vq_class .* in_africa;

%% (7.2) Plot the map showing classification 
figure(211)
pcolor(xq,yq, Vq_class)
shading flat
% c1 = colorbar;
% c1.Title.String = 'index';
colormap([0,0.4471,0.7412;0.4667,0.6745,0.1882;0.8510,0.3255,0.0980;0.9294,0.6941,0.1255])

hold on
% plot continent
plot(coastline, 'EdgeColor', 'k', 'LineWidth', 4, 'FaceAlpha', 0);  % plot the African boundary
plot(madgline, 'EdgeColor', 'k', 'LineWidth', 4, 'FaceAlpha', 0);  % plot the African boundary
% plot African cratons
plot(AC.Congo(:,2), AC.Congo(:,1),':', 'color', [0.85 0.85 0.85], 'lineWidth', 2);
plot(AC.Kala(:,2), AC.Kala(:,1),':', 'color', [0.85 0.85 0.85], 'lineWidth', 2);
plot(AC.Sahara(:,2), AC.Sahara(:,1),':', 'color', [0.85 0.85 0.85], 'lineWidth', 2);
plot(AC.Tanza(:,2), AC.Tanza(:,1),':', 'color', [0.85 0.85 0.85], 'lineWidth', 2);
plot(AC.West(:,2), AC.West(:,1), ':','color', [0.85 0.85 0.85], 'lineWidth', 2);
plot(AC.Zim(:,2), AC.Zim(:,1), ':','color', [0.85 0.85 0.85], 'lineWidth', 2);
% colorbar

% plot crustal regions
plot(GEO.MobileBelts.Other.Cap_belt(:,2), GEO.MobileBelts.Other.Cap_belt(:,1), 'k','lineWidth', 1)
plot(GEO.MobileBelts.Other.Kibaran_belt(:,2), GEO.MobileBelts.Other.Kibaran_belt(:,1), 'k','lineWidth', 1)
plot(GEO.MobileBelts.Other.Namaquanatal(:,2), GEO.MobileBelts.Other.Namaquanatal(:,1), 'k','lineWidth', 1)
plot(GEO.MobileBelts.Other.Ruwenzky(:,2), GEO.MobileBelts.Other.Ruwenzky(:,1), 'k','lineWidth', 1)
plot(GEO.MobileBelts.Other.Mauritania(:,2), GEO.MobileBelts.Other.Mauritania(:,1), 'k','lineWidth', 1)
plot(GEO.MobileBelts.WestAfrica_belt(:,2), GEO.MobileBelts.WestAfrica_belt(:,1), 'k','lineWidth', 1)
plot(GEO.MobileBelts.Ouban_Damara.Damara(:,2), GEO.MobileBelts.Ouban_Damara.Damara(:,1), 'k','lineWidth', 1)
plot(GEO.MobileBelts.Ouban_Damara.Oubangides(:,2), GEO.MobileBelts.Ouban_Damara.Oubangides(:,1), 'k','lineWidth', 1)

plot(GEO.Orogens.AtlasMountains(:,2), GEO.Orogens.AtlasMountains(:,1), 'k','lineWidth', 1)
plot(GEO.Orogens.EastAfrica.EastAfrica(:,2), GEO.Orogens.EastAfrica.EastAfrica(:,1), 'k','lineWidth', 1)
plot(GEO.Orogens.EastAfrica.NorthMozamb(:,2), GEO.Orogens.EastAfrica.NorthMozamb(:,1), 'k','lineWidth', 1)

plot(GEO.Basins.Congo(:,2), GEO.Basins.Congo(:,1), 'k','lineWidth', 1)
plot(GEO.Basins.Tindouf.Tindouf(:,2), GEO.Basins.Tindouf.Tindouf(:,1), 'k','lineWidth', 1)
plot(GEO.Basins.Tindouf.MauriSenba(:,2), GEO.Basins.Tindouf.MauriSenba(:,1), 'k','lineWidth', 1)
plot(GEO.Basins.Taoudeni.B1(:,2), GEO.Basins.Taoudeni.B1(:,1), 'k','lineWidth', 1)

plot(GEO.Archean.WestAfrica.LeoManShield(:,2), GEO.Archean.WestAfrica.LeoManShield(:,1), 'k','lineWidth', 1)
plot(GEO.Archean.WestAfrica.Reguibatshield(:,2), GEO.Archean.WestAfrica.Reguibatshield(:,1), 'k','lineWidth', 1)
plot(GEO.Archean.Tanzania(:,2), GEO.Archean.Tanzania(:,1), 'k','lineWidth', 1)
plot(GEO.Archean.CongoAll.AC(:,2), GEO.Archean.CongoAll.AC(:,1), 'k','lineWidth', 1)
plot(GEO.Archean.CongoAll.Bomu(:,2), GEO.Archean.CongoAll.Bomu(:,1), 'k','lineWidth', 1)
% plot(GEO.Archean.CongoAll.Congo(:,2), GEO.Archean.CongoAll.Congo(:,1), 'k','lineWidth', 1)
plot(GEO.Archean.CongoAll.Kasai(:,2), GEO.Archean.CongoAll.Kasai(:,1), 'k','lineWidth', 1)
plot(GEO.Archean.CongoAll.Ntem(:,2), GEO.Archean.CongoAll.Ntem(:,1), 'k','lineWidth', 1)
plot(GEO.Archean.SaharaMeta.A1(:,2), GEO.Archean.SaharaMeta.A1(:,1), 'k','lineWidth', 1)
plot(GEO.Archean.SaharaMeta.A2(:,2), GEO.Archean.SaharaMeta.A2(:,1), 'k','lineWidth', 1)
plot(GEO.Archean.SaharaMeta.A3(:,2), GEO.Archean.SaharaMeta.A3(:,1), 'k','lineWidth', 1)
plot(GEO.Archean.SaharaMeta.A4(:,2), GEO.Archean.SaharaMeta.A4(:,1), 'k','lineWidth', 1)
plot(GEO.Archean.WestAfrica_MZ.Dahomeyshield(:,2), GEO.Archean.WestAfrica_MZ.Dahomeyshield(:,1), 'k','lineWidth', 1)
plot(GEO.Archean.WestAfrica_MZ.Tuaregshield(:,2), GEO.Archean.WestAfrica_MZ.Tuaregshield(:,1), 'k','lineWidth', 1)
plot(GEO.Archean.KaapVaal(:,2), GEO.Archean.KaapVaal(:,1), 'k','lineWidth', 1)

% xlabel('Longitude')
% ylabel('Latitude')
minlon = -20;
maxlon = 55;
minlat = -36;
maxlat = 40;
xlim([minlon maxlon]);
ylim([minlat maxlat]);
axesm('mercator','maplatlimit',[minlat maxlat],'maplonlimit',[minlon maxlon]);
xticks(10*ceil(minlon/10):10:10*floor(maxlon/10));
yticks(10*ceil(minlat/10):10:10*floor(maxlat/10));
% title('K-Mean Clustering')
set(gca,'fontsize', 14) 
set(gcf, 'color', 'w')


%% (7.2-2) manually create a legend for class types
figure(221)
bar(0,1,'FaceColor', [0,0.4471,0.7412])
hold on
bar(0,1,'FaceColor', [0.9294,0.6941,0.1255])
bar(0,1,'FaceColor', [0.4667,0.6745,0.1882])
bar(0,1,'FaceColor', [0.8510,0.3255,0.0980])

hleg1 = legend({'C1','C2','C3','C4'}, 'fontsize', 12);
set(hleg1,'position',[0.2 0.2 0.1 0.18])

%% (7.3) plot the vel_depth distribution
%%%%%%%%% Note: the the class number not necessarly equal to class idx!!!!!!!!!!
t = tiledlayout(4,1,'TileSpacing','tight','Padding','tight');
nexttile
C1_vel = vel1D_ADAMA(Afrclass == 1, :)/1000;
C1_avg = mean(C1_vel(:, 1:81));
hold on
for i = 1:size(C1_vel, 1)
    plot(C1_vel(i, 1:81),0:0.5:40, 'Color', [0.651, 0.651, 0.651])
end

plot(C1_avg,0:0.5:40, 'Color', [0,0.4471,0.7412], 'LineWidth', 5)
hold off
text(2.1,35,['\sigma^2=', num2str(round(sum(var(C1_vel(:,1:81))),2))],'FontSize',14)
text(2.1,25,[num2str(round(size(C1_vel, 1)/2381*100, 1)), '%'],'FontSize',14)
xlim([2, 5])
xticklabels({})
set(gca, 'YDir','reverse','fontsize', 14)

nexttile
C2_vel = vel1D_ADAMA(Afrclass == 4, :)/1000;
C2_avg = mean(C2_vel(:, 1:81));
hold on
for i = 1:size(C2_vel, 1)
    plot(C2_vel(i, 1:81),0:0.5:40, 'Color', [0.651, 0.651, 0.651])
end
plot(C2_avg,0:0.5:40, 'Color', [0.9294,0.6941,0.1255], 'LineWidth', 5)
hold off
text(2.1,35,['\sigma^2=', num2str(round(sum(var(C2_vel(:,1:81))),2))],'FontSize',14)
text(2.1,25,[num2str(round(size(C2_vel, 1)/2381*100, 1)), '%'],'FontSize',14)
xlim([2, 5])
xticklabels({})
set(gca, 'YDir','reverse','fontsize', 14)

nexttile
C3_vel = vel1D_ADAMA(Afrclass == 2, :)/1000;
C3_avg = mean(C3_vel(:, 1:81));
hold on
for i = 1:size(C3_vel, 1)
    plot(C3_vel(i, 1:81),0:0.5:40, 'Color', [0.651, 0.651, 0.651])
end
plot(C3_avg,0:0.5:40, 'Color', [0.4667,0.6745,0.1882], 'LineWidth', 5)
hold off
text(2.1,35,['\sigma^2=', num2str(round(sum(var(C3_vel(:,1:81))),2))],'FontSize',14)
text(2.1,25,[num2str(round(size(C3_vel, 1)/2381*100, 1)), '%'],'FontSize',14)
xlim([2, 5])
xticklabels({})
set(gca, 'YDir','reverse','fontsize', 14)

nexttile
C4_vel = vel1D_ADAMA(Afrclass == 3, :)/1000;
C4_avg = mean(C4_vel(:, 1:81));
hold on
for i = 1:size(C4_vel, 1)
    plot(C4_vel(i, 1:81),0:0.5:40, 'Color', [0.651, 0.651, 0.651])
end
plot(C4_avg,0:0.5:40, 'Color', [0.8510,0.3255,0.0980], 'LineWidth', 5)
hold off
text(2.1,35,['\sigma^2=', num2str(round(sum(var(C4_vel(:,1:81))),2))],'FontSize',14)
text(2.1,25,[num2str(round(size(C4_vel, 1)/2381*100, 1)), '%'],'FontSize',14)
set(gca, 'YDir','reverse','fontsize', 14)
xlim([2, 5])
xlabel('Velocity (km/s)')
ylabel('Depth (km)')


set(gcf, 'color', 'w','position',[200, 200, 200, 500])

%% (7.4) plot the vel_depth distribution subclass (k=2 per each vel class)

Subclass = readmatrix('/Users/sxue3/Documents/BayMap_Figures/Data/Sequencer/AfricaClass_k4_subk2.csv');

t = tiledlayout(4,1,'TileSpacing','tight','Padding','tight');
nexttile
C11_vel = vel1D_ADAMA(Subclass == 1, :)/1000;
hold on
for i = 1:size(C11_vel, 1)
    plot(C11_vel(i, 1:81),0:0.5:40, 'Color', [0.651, 0.651, 0.651])
end


C12_vel = vel1D_ADAMA(Subclass == 2, :)/1000;
for i = 1:size(C12_vel, 1)
    plot(C12_vel(i, 1:81),0:0.5:40, 'Color', [0.8 0.8 0.8])
end
plot(mean(C11_vel(:, 1:81)),0:0.5:40, 'Color', [0,0.4471,0.7412], 'LineWidth', 5)
plot(mean(C12_vel(:, 1:81)),0:0.5:40, 'Color', [0.302, 0.745, 0.933], 'LineWidth', 5)
hold off
% text(2.1,35,['\sigma^2=', num2str(round(sum(var(C1_vel(:,1:81))),2))],'FontSize',14)
% text(2.1,25,[num2str(round(size(C1_vel, 1)/2381*100, 1)), '%'],'FontSize',14)
xlim([2, 5])
set(gca, 'YDir','reverse','fontsize', 14)

nexttile
C21_vel = vel1D_ADAMA(Subclass == 3, :)/1000;

hold on
for i = 1:size(C21_vel, 1)
    plot(C21_vel(i, 1:81),0:0.5:40, 'Color', [0.651, 0.651, 0.651])
end


C22_vel = vel1D_ADAMA(Subclass == 4, :)/1000;
for i = 1:size(C22_vel, 1)
    plot(C22_vel(i, 1:81),0:0.5:40, 'Color', [0.8 0.8 0.8])
end
plot(mean(C21_vel(:, 1:81)),0:0.5:40, 'Color', [0.9294,0.6941,0.1255], 'LineWidth', 5)
plot(mean(C22_vel(:, 1:81)),0:0.5:40, 'Color', [0.98, 0.863, 0.392], 'LineWidth', 5)
hold off
% text(2.1,35,['\sigma^2=', num2str(round(sum(var(C2_vel(:,1:81))),2))],'FontSize',14)
% text(2.1,25,[num2str(round(size(C2_vel, 1)/2381*100, 1)), '%'],'FontSize',14)
xlim([2, 5])
set(gca, 'YDir','reverse','fontsize', 14)

nexttile
C31_vel = vel1D_ADAMA(Subclass == 5, :)/1000;
hold on
for i = 1:size(C31_vel, 1)
    plot(C31_vel(i, 1:81),0:0.5:40, 'Color', [0.651, 0.651, 0.651])
end


C32_vel = vel1D_ADAMA(Subclass == 6, :)/1000;
hold on
for i = 1:size(C32_vel, 1)
    plot(C32_vel(i, 1:81),0:0.5:40, 'Color', [0.8 0.8 0.8])
end
plot(mean(C31_vel(:, 1:81)),0:0.5:40, 'Color', [0.4667,0.6745,0.1882], 'LineWidth', 5)
plot(mean(C32_vel(:, 1:81)),0:0.5:40, 'Color', [0.714, 0.859, 0.129], 'LineWidth', 5)
hold off
% text(2.1,35,['\sigma^2=', num2str(round(sum(var(C3_vel(:,1:81))),2))],'FontSize',14)
% text(2.1,25,[num2str(round(size(C3_vel, 1)/2381*100, 1)), '%'],'FontSize',14)
xlim([2, 5])
set(gca, 'YDir','reverse','fontsize', 14)

nexttile
C41_vel = vel1D_ADAMA(Subclass == 7, :)/1000;
hold on
for i = 1:size(C41_vel, 1)
    plot(C41_vel(i, 1:81),0:0.5:40, 'Color', [0.651, 0.651, 0.651])
end


C42_vel = vel1D_ADAMA(Subclass == 8, :)/1000;
hold on
for i = 1:size(C42_vel, 1)
    plot(C42_vel(i, 1:81),0:0.5:40, 'Color', [0.8 0.8 0.8])
end
plot(mean(C41_vel(:, 1:81)),0:0.5:40, 'Color', [0.8510,0.3255,0.0980], 'LineWidth', 5)
plot(mean(C42_vel(:, 1:81)),0:0.5:40, 'Color', [0.989, 0.467, 0.118], 'LineWidth', 5)
hold off
% text(2.1,35,['\sigma^2=', num2str(round(sum(var(C4_vel(:,1:81))),2))],'FontSize',14)
% text(2.1,25,[num2str(round(size(C4_vel, 1)/2381*100, 1)), '%'],'FontSize',14)
xlim([2, 5])
set(gca, 'YDir','reverse','fontsize', 14)
% xlabel('Velocity (km/s)')
xlabel('Velocity Profiles')
ylabel('Depth (km)')

set(gcf, 'color', 'w','position',[200, 200, 200, 500])

%% (7.5) Plot Africa subclass on the map

figure(742)
hold on
for i = 1:2381
    if Subclass(i) == 1
        plot(dslon(i), dslat(i), '.', 'MarkerSize',20, 'Color',[0,0.4471,0.7412])
    elseif Subclass(i) == 2
        plot(dslon(i), dslat(i), '.', 'MarkerSize',20, 'Color',[0.302, 0.745, 0.933])
    elseif Subclass(i) == 3
        plot(dslon(i), dslat(i), '.', 'MarkerSize',20, 'Color',[0.9294,0.6941,0.1255])
    elseif Subclass(i) == 4
        plot(dslon(i), dslat(i), '.', 'MarkerSize',20, 'Color',[0.98, 0.863, 0.392])
    elseif Subclass(i) == 5
        plot(dslon(i), dslat(i), '.', 'MarkerSize',20, 'Color',[0.4667,0.6745,0.1882])
    elseif Subclass(i) == 6
        plot(dslon(i), dslat(i), '.', 'MarkerSize',20, 'Color',[0.714, 0.859, 0.129])
    elseif Subclass(i) == 7
        plot(dslon(i), dslat(i), '.', 'MarkerSize',20, 'Color',[0.8510,0.3255,0.0980])
    else
        plot(dslon(i), dslat(i), '.', 'MarkerSize',20, 'Color',[0.989, 0.467, 0.118])
    end

end

% plot continent
plot(coastline, 'EdgeColor', 'k', 'LineWidth', 4, 'FaceAlpha', 0);  % plot the African boundary
plot(madgline, 'EdgeColor', 'k', 'LineWidth', 4, 'FaceAlpha', 0);  % plot the African boundary

% plot crustal regions
plot(GEO.MobileBelts.Other.Cap_belt(:,2), GEO.MobileBelts.Other.Cap_belt(:,1), 'k','lineWidth', 1)
plot(GEO.MobileBelts.Other.Kibaran_belt(:,2), GEO.MobileBelts.Other.Kibaran_belt(:,1), 'k','lineWidth', 1)
plot(GEO.MobileBelts.Other.Namaquanatal(:,2), GEO.MobileBelts.Other.Namaquanatal(:,1), 'k','lineWidth', 1)
plot(GEO.MobileBelts.Other.Ruwenzky(:,2), GEO.MobileBelts.Other.Ruwenzky(:,1), 'k','lineWidth', 1)
plot(GEO.MobileBelts.Other.Mauritania(:,2), GEO.MobileBelts.Other.Mauritania(:,1), 'k','lineWidth', 1)
plot(GEO.MobileBelts.WestAfrica_belt(:,2), GEO.MobileBelts.WestAfrica_belt(:,1), 'k','lineWidth', 1)
plot(GEO.MobileBelts.Ouban_Damara.Damara(:,2), GEO.MobileBelts.Ouban_Damara.Damara(:,1), 'k','lineWidth', 1)
plot(GEO.MobileBelts.Ouban_Damara.Oubangides(:,2), GEO.MobileBelts.Ouban_Damara.Oubangides(:,1), 'k','lineWidth', 1)

plot(GEO.Orogens.AtlasMountains(:,2), GEO.Orogens.AtlasMountains(:,1), 'k','lineWidth', 1)
plot(GEO.Orogens.EastAfrica.EastAfrica(:,2), GEO.Orogens.EastAfrica.EastAfrica(:,1), 'k','lineWidth', 1)
plot(GEO.Orogens.EastAfrica.NorthMozamb(:,2), GEO.Orogens.EastAfrica.NorthMozamb(:,1), 'k','lineWidth', 1)

plot(GEO.Basins.Congo(:,2), GEO.Basins.Congo(:,1), 'k','lineWidth', 1)
plot(GEO.Basins.Tindouf.Tindouf(:,2), GEO.Basins.Tindouf.Tindouf(:,1), 'k','lineWidth', 1)
plot(GEO.Basins.Tindouf.MauriSenba(:,2), GEO.Basins.Tindouf.MauriSenba(:,1), 'k','lineWidth', 1)
plot(GEO.Basins.Taoudeni.B1(:,2), GEO.Basins.Taoudeni.B1(:,1), 'k','lineWidth', 1)

plot(GEO.Archean.WestAfrica.LeoManShield(:,2), GEO.Archean.WestAfrica.LeoManShield(:,1), 'k','lineWidth', 1)
plot(GEO.Archean.WestAfrica.Reguibatshield(:,2), GEO.Archean.WestAfrica.Reguibatshield(:,1), 'k','lineWidth', 1)
plot(GEO.Archean.Tanzania(:,2), GEO.Archean.Tanzania(:,1), 'k','lineWidth', 1)
plot(GEO.Archean.CongoAll.AC(:,2), GEO.Archean.CongoAll.AC(:,1), 'k','lineWidth', 1)
plot(GEO.Archean.CongoAll.Bomu(:,2), GEO.Archean.CongoAll.Bomu(:,1), 'k','lineWidth', 1)
% plot(GEO.Archean.CongoAll.Congo(:,2), GEO.Archean.CongoAll.Congo(:,1), 'k','lineWidth', 1)
plot(GEO.Archean.CongoAll.Kasai(:,2), GEO.Archean.CongoAll.Kasai(:,1), 'k','lineWidth', 1)
plot(GEO.Archean.CongoAll.Ntem(:,2), GEO.Archean.CongoAll.Ntem(:,1), 'k','lineWidth', 1)
plot(GEO.Archean.SaharaMeta.A1(:,2), GEO.Archean.SaharaMeta.A1(:,1), 'k','lineWidth', 1)
plot(GEO.Archean.SaharaMeta.A2(:,2), GEO.Archean.SaharaMeta.A2(:,1), 'k','lineWidth', 1)
plot(GEO.Archean.SaharaMeta.A3(:,2), GEO.Archean.SaharaMeta.A3(:,1), 'k','lineWidth', 1)
plot(GEO.Archean.SaharaMeta.A4(:,2), GEO.Archean.SaharaMeta.A4(:,1), 'k','lineWidth', 1)
plot(GEO.Archean.WestAfrica_MZ.Dahomeyshield(:,2), GEO.Archean.WestAfrica_MZ.Dahomeyshield(:,1), 'k','lineWidth', 1)
plot(GEO.Archean.WestAfrica_MZ.Tuaregshield(:,2), GEO.Archean.WestAfrica_MZ.Tuaregshield(:,1), 'k','lineWidth', 1)
plot(GEO.Archean.KaapVaal(:,2), GEO.Archean.KaapVaal(:,1), 'k','lineWidth', 1)

% xlabel('Longitude')
% ylabel('Latitude')
minlon = -20;
maxlon = 55;
minlat = -36;
maxlat = 40;
xlim([minlon maxlon]);
ylim([minlat maxlat]);
axesm('mercator','maplatlimit',[minlat maxlat],'maplonlimit',[minlon maxlon]);
xticks(10*ceil(minlon/10):10:10*floor(maxlon/10));
yticks(10*ceil(minlat/10):10:10*floor(maxlat/10));
% title('K-Mean Clustering')
set(gca,'fontsize', 14) 
set(gcf, 'color', 'w')


%% (7.5-2) manually create a legend for subclass types
figure(221)
hold on
bar(0,1,'FaceColor', [0,0.4471,0.7412])
bar(0,1,'FaceColor', [0.302, 0.745, 0.933])
bar(0,1,'FaceColor', [0.9294,0.6941,0.1255])
bar(0,1,'FaceColor', [0.98, 0.863, 0.392])
bar(0,1,'FaceColor', [0.4667,0.6745,0.1882])
bar(0,1,'FaceColor', [0.714, 0.859, 0.129])
bar(0,1,'FaceColor', [0.8510,0.3255,0.0980])
bar(0,1,'FaceColor', [0.989, 0.467, 0.118])

hleg1 = legend({'C1-1','C1-2','C2-1','C2-2','C3-1','C3-2','C4-1','C4-2'}, 'fontsize', 12);
set(hleg1,'position',[0.2 0.2 0.1 0.18])


%%
temp = [C1_vel; C2_vel; C3_vel; C4_vel];
pcolor(temp)
shading flat

%% (7.6) compare avg vel of each class

plot(mean(C1_vel(:, 1:81)),0:0.5:40, 'Color', [0,0.4471,0.7412], 'LineWidth', 5)
hold on
plot(mean(C2_vel(:, 1:81)),0:0.5:40, 'Color', [0.9294,0.6941,0.1255], 'LineWidth', 5)
plot(mean(C3_vel(:, 1:81)),0:0.5:40, 'Color', [0.4667,0.6745,0.1882], 'LineWidth', 5)
plot(mean(C4_vel(:, 1:81)),0:0.5:40, 'Color', [0.8510,0.3255,0.0980], 'LineWidth', 5)
set(gca, 'YDir','reverse','fontsize', 14)

%% (9.1) Plot the classification using correlation cutoff


CUTs = readmatrix('/Users/sxue3/Documents/BayMap_Figures/Data/Sequencer/boundary_locations_05.txt');
sqe_idx = IDX(:,2);
corr_class = ones(size(sqe_idx)) * (length(CUTs) + 1);

curr_class = 1;
for i = 1:length(corr_class)
    corr_class(i) = curr_class;

    if curr_class == length(CUTs)
        break
    end

    if sqe_idx(i) == CUTs(curr_class)
        curr_class = curr_class + 1;
    end
end

%% interpolate the map
xgrid = -20:0.1:55;    % lon
ygrid = -35:0.1:40;    % lat
[xq, yq] = meshgrid(xgrid, ygrid);
F = scatteredInterpolant(dslon(sqe_idx+1),dslat(sqe_idx+1),corr_class, 'nearest', 'none');
Vq_class = F(xq,yq);

Vq_class = Vq_class .* in_africa;


%% (9.2) Plot the map showing index 
pcolor(xq,yq, Vq_class)
shading flat
c1 = colorbar;
c1.Title.String = 'Class';
colormap('jet')

hold on
% plot continent
plot(coastline, 'EdgeColor', 'k', 'LineWidth', 4, 'FaceAlpha', 0);  % plot the African boundary
plot(madgline, 'EdgeColor', 'k', 'LineWidth', 4, 'FaceAlpha', 0);  % plot the African boundary
% plot African cratons
plot(AC.Congo(:,2), AC.Congo(:,1), 'color', 'k', 'lineWidth', 2);
plot(AC.Kala(:,2), AC.Kala(:,1), 'color', 'k', 'lineWidth', 2);
plot(AC.Sahara(:,2), AC.Sahara(:,1), 'color', 'k', 'lineWidth', 2);
plot(AC.Tanza(:,2), AC.Tanza(:,1), 'color', 'k', 'lineWidth', 2);
plot(AC.West(:,2), AC.West(:,1), 'color', 'k', 'lineWidth', 2);
plot(AC.Zim(:,2), AC.Zim(:,1), 'color', 'k', 'lineWidth', 2);

xlabel('Longitude')
ylabel('Latitude')
minlon = -20;
maxlon = 55;
minlat = -36;
maxlat = 40;
xlim([minlon maxlon]);
ylim([minlat maxlat]);
axesm('mercator','maplatlimit',[minlat maxlat],'maplonlimit',[minlon maxlon]);
xticks(10*ceil(minlon/10):10:10*floor(maxlon/10));
yticks(10*ceil(minlat/10):10:10*floor(maxlat/10));
title('Sequencer Index')
set(gca,'fontsize', 14) 
set(gcf, 'color', 'w')
