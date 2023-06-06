%%%%% Fig 4: average velocity maps, standard deviation maps, skewness maps, and Litho maps
%%%%% Siyu Xue -- Mar 10. 2023
addpath('/Users/sxue3/Documents/BayMap_Figures/Utils')
addpath('/Users/sxue3/Documents/ADAMA_Figures')

%% Read in data
clear

% load tha map values
load('/Users/sxue3/Documents/BayMap_Figures/Data/LoveAvgMap/L35_P_maps.mat');

% don't change the following two paths (plotpara.mat is the same for all maps)
longrid = load('/Users/sxue3/Documents/BayMap_Figures/Data/RayleighPData/plotpara.mat').longrid;
latgrid = load('/Users/sxue3/Documents/BayMap_Figures/Data/RayleighPData/plotpara.mat').latgrid;

minlat = min(latgrid(:,1)); maxlat = max(latgrid(:,1));
minlon = min(longrid(1,:)); maxlon = max(longrid(1,:));

% Read in the percentage error 
err_path = '/Users/sxue3/Documents/BayMap_Figures/Data/map_error/observation_errors_L35P.txt';
A = readtable(err_path,'FileType','text');
percenErr= A{:,'Var8'};
clear A

% load African country boundaries
S1 = load('/Users/sxue3/Documents/BayMap_Figures/Data/AfrCountry.mat').S1;
madgline = polyshape((S1(448).X).',(S1(448).Y).');  % load Madgascar's Big Island

% load the coast
afcoast = load('/Users/sxue3/Documents/BayMap_Figures/Data/AfricaCoast.mat');
aflon = wrapTo180(afcoast.XY(:, 1));
aflat = afcoast.XY(:,2);
coastline = polyshape(aflon, aflat);

% load the Tectonic Plate Boundaries
Tplates = shaperead('./Data/Tectonic_PB.shp');

% load cratons
load('/Users/sxue3/Documents/BayMap_Figures/Data/AfricanCratons.mat');

% Load all points generated using equal distance (for Litho)
load('/Users/sxue3/Documents/BayMap_Figures/Data/AfrDownSampleDots.mat');

% Load the Litho velocity
load('/Users/sxue3/Documents/BayMap_Figures/Data/Litho_dwsample.mat')

Periods = [40 35 30 25 20 15 12 10 8 6 5];
dslove = LithoLove(:,2)./1000;  % Don't forget the change index for different period

set(0,'defaultfigurecolor',[1 1 1])

%% Create a mask to seperate region within Africa
afrmask = zeros(297, 324);
mapsize = size(afrmask);
for xi = 1:mapsize(2)
    loni = longrid(1,xi);
    for yi = 1:mapsize(1)
        lati = latgrid(yi,1);

        if isinterior(coastline, loni, lati) == 1 || isinterior(madgline, loni, lati) == 1
            afrmask(yi, xi) = 1;
        end
    end
end

afrmask = logical(afrmask - 1);  % points not in Africa == 1

%% (a1) Plot the average velocity map
ax1 = subplot(2,3,[1,2,4,5]);
avgmap(afrmask) = nan;
pcolor(ax1, longrid,latgrid,avgmap);
shading flat;
caxis([mean(avgmap, [1,2], 'omitnan') - 2*std(avgmap, 0, [1,2], 'omitnan'), ...
    mean(avgmap, [1,2], 'omitnan') + 2*std(avgmap, 0, [1,2], 'omitnan')]); 
cb1 = colorbar;
title(cb1, 'Velocity(km/s)', 'FontSize', 12)

% customize the colormap colors
RedMap = [linspace(0.98, 1, 128)', linspace(0.16, 1, 128)', linspace(0.16, 1, 128)'];
BlueMap = [linspace(1, 0.16, 128)', linspace(1, 0.57, 128)', linspace(1, 0.98, 128)'];
RWBMap = [RedMap; BlueMap];
colormap(ax1, RWBMap);

hold on

% Plot the African boundary
plot(coastline, 'EdgeColor', 'k', 'LineWidth', 3.5, 'FaceAlpha', 0);  
% Plot Madgascar's Big Island
plot(madgline.Vertices(1:50:end, 1), madgline.Vertices(1:50:end, 2),'k', 'LineWidth', 3.5);
% Plot tectonic plate outlines
for ip = 1:117
    plot(Tplates(ip).X, Tplates(ip).Y, 'LineWidth', 3, 'Color',[0.4 0.4 0.4])
end

% Plot African cratons
plot(AC.Congo(:,2), AC.Congo(:,1), 'color', [0.6 0.6 0.6], 'lineWidth', 2);
plot(AC.Kala(:,2), AC.Kala(:,1), 'color', [0.6 0.6 0.6], 'lineWidth', 2);
plot(AC.Sahara(:,2), AC.Sahara(:,1), 'color', [0.6 0.6 0.6], 'lineWidth', 2);
plot(AC.Tanza(:,2), AC.Tanza(:,1), 'color', [0.6 0.6 0.6], 'lineWidth', 2);
plot(AC.West(:,2), AC.West(:,1), 'color', [0.6 0.6 0.6], 'lineWidth', 2);
plot(AC.Zim(:,2), AC.Zim(:,1), 'color', [0.6 0.6 0.6], 'lineWidth', 2);

xlim([minlon maxlon]);
ylim([minlat maxlat]);
axesm('mercator','maplatlimit',[minlat maxlat],'maplonlimit',[minlon maxlon]);
xticks(10*ceil(minlon/10):10:10*floor(maxlon/10));
yticks(10*ceil(minlat/10):10:10*floor(maxlat/10));

annotation('textbox',[0.133, 0.875, 0.06, 0.05], 'String', ...
    'L35-P','EdgeColor','k', 'BackgroundColor','w', 'FontSize', 18)
hold off


% (a2) Plot the standard deviation map
ax2 = axes('position',[0.6, 0.54, 0.38, 0.38]);  % if plot skewness not Litho
% ax2 = axes('position',[0.6, 0.11, 0.38, 0.38]);
sdevmap(afrmask) = nan;
pcolor(ax2, longrid,latgrid,sdevmap);
shading flat; 
cb2 = colorbar;
%set(cb2,'position',[0.92, 0.18, 0.02, 0.5])
title(cb2, 'SD', 'FontSize', 12)
colormap(ax2, flipud(hot));
caxis([0 1]); 

hold on

% Plot the African boundary
plot(coastline, 'EdgeColor', 'k', 'LineWidth', 2.5, 'FaceAlpha', 0);
% Plot Madgascar's Big Island
plot(madgline.Vertices(1:50:end, 1), madgline.Vertices(1:50:end, 2),'k', 'LineWidth', 2.5);
% Plot tectonic plate outlines
for ip = 1:117
    plot(Tplates(ip).X, Tplates(ip).Y, 'LineWidth', 2, 'Color',[0.4 0.4 0.4])
end

% Plot African cratons
plot(AC.Congo(:,2), AC.Congo(:,1), 'color', [0.7 0.7 0.7], 'lineWidth', 2);
plot(AC.Kala(:,2), AC.Kala(:,1), 'color', [0.7 0.7 0.7], 'lineWidth', 2);
plot(AC.Sahara(:,2), AC.Sahara(:,1), 'color', [0.7 0.7 0.7], 'lineWidth', 2);
plot(AC.Tanza(:,2), AC.Tanza(:,1), 'color', [0.7 0.7 0.7], 'lineWidth', 2);
plot(AC.West(:,2), AC.West(:,1), 'color', [0.7 0.7 0.7], 'lineWidth', 2);
plot(AC.Zim(:,2), AC.Zim(:,1), 'color', [0.7 0.7 0.7], 'lineWidth', 2);

xlim([minlon maxlon]);
ylim([minlat maxlat]);
axesm('mercator','maplatlimit',[minlat maxlat],'maplonlimit',[minlon maxlon]);

hold off

% % (a3) Plot litho velocity map
% ax3 = axes('position',[0.058, 0.11, 0.33, 0.33]);
% dsvq = griddata(dslon,dslat,dslove,longrid,latgrid, 'nearest');
% dsvq(afrmask) = nan;
% pcolor(ax3,longrid,latgrid, dsvq)
% shading flat;
% 
% caxis([mean(avgmap, [1,2], 'omitnan') - 1*std(avgmap, 0, [1,2], 'omitnan'), ...
%     mean(avgmap, [1,2], 'omitnan') + 1*std(avgmap, 0, [1,2], 'omitnan')]); 
% colormap(ax3, RWBMap);
% 
% hold on
% 
% % Plot the African boundary
% plot(coastline, 'EdgeColor', 'k', 'LineWidth', 2, 'FaceAlpha', 0);
% % Plot Madgascar's Big Island
% plot(madgline.Vertices(1:50:end, 1), madgline.Vertices(1:50:end, 2),'k', 'LineWidth', 2);
% % Plot tectonic plate outlines
% for ip = 1:117
%     plot(Tplates(ip).X, Tplates(ip).Y, 'LineWidth', 2, 'Color',[0.4 0.4 0.4])
% end
% 
% % Plot African cratons
% plot(AC.Congo(:,2), AC.Congo(:,1), 'color', [0.6 0.6 0.6], 'lineWidth', 1.5);
% plot(AC.Kala(:,2), AC.Kala(:,1), 'color', [0.6 0.6 0.6], 'lineWidth', 1.5);
% plot(AC.Sahara(:,2), AC.Sahara(:,1), 'color', [0.6 0.6 0.6], 'lineWidth', 1.5);
% plot(AC.Tanza(:,2), AC.Tanza(:,1), 'color', [0.6 0.6 0.6], 'lineWidth', 1.5);
% plot(AC.West(:,2), AC.West(:,1), 'color', [0.6 0.6 0.6], 'lineWidth', 1.5);
% plot(AC.Zim(:,2), AC.Zim(:,1), 'color', [0.6 0.6 0.6], 'lineWidth', 1.5);
% 
% xlim([minlon maxlon]);
% ylim([minlat maxlat]);
% axesm('mercator','maplatlimit',[minlat maxlat],'maplonlimit',[minlon maxlon]);
% 
% hold off

% (a4) show the % of observations have error greater than 20%
badtotal = sum(percenErr > 20)/length(percenErr) * 100;
badtotal = sprintf("%.2f",badtotal);
badtotal = convertStringsToChars(badtotal);

annotation('textbox',[0.133, 0.05, 0.5, 0.05], 'String', ...
    [badtotal, '% of the total observations exceed 20 percentage error'],...
    'EdgeColor','k', 'BackgroundColor','w', 'FontSize', 18)


% % (a5) Plot the (ADAMA - Litho) velocity map
% ax5 = axes('position',[0.6, 0.54, 0.38, 0.38]);
% diffmap = abs(avgmap - dsvq);
% pcolor(ax5,longrid,latgrid, diffmap)
% shading flat;
% 
% cb5 = colorbar;
% title(cb5, 'Velocity(km/s)', 'FontSize', 12)
% colormap(ax5, flipud(hot));
% caxis([0, 1])
% 
% hold on
% 
% % Plot the African boundary
% plot(coastline, 'EdgeColor', 'k', 'LineWidth', 2.5, 'FaceAlpha', 0);
% % Plot Madgascar's Big Island
% plot(madgline.Vertices(1:50:end, 1), madgline.Vertices(1:50:end, 2),'k', 'LineWidth', 2.5);
% % Plot tectonic plate outlines
% for ip = 1:117
%     plot(Tplates(ip).X, Tplates(ip).Y, 'LineWidth', 2, 'Color',[0.4 0.4 0.4])
% end
% 
% % Plot African cratons
% plot(AC.Congo(:,2), AC.Congo(:,1), 'color', [0.7 0.7 0.7], 'lineWidth', 2);
% plot(AC.Kala(:,2), AC.Kala(:,1), 'color', [0.7 0.7 0.7], 'lineWidth', 2);
% plot(AC.Sahara(:,2), AC.Sahara(:,1), 'color', [0.7 0.7 0.7], 'lineWidth', 2);
% plot(AC.Tanza(:,2), AC.Tanza(:,1), 'color', [0.7 0.7 0.7], 'lineWidth', 2);
% plot(AC.West(:,2), AC.West(:,1), 'color', [0.7 0.7 0.7], 'lineWidth', 2);
% plot(AC.Zim(:,2), AC.Zim(:,1), 'color', [0.7 0.7 0.7], 'lineWidth', 2);
% 
% xlim([minlon maxlon]);
% ylim([minlat maxlat]);
% axesm('mercator','maplatlimit',[minlat maxlat],'maplonlimit',[minlon maxlon]);
% 
% hold off

% (a6) Plot the distribution of skewness
ax4 = axes('position',[0.6, 0.11, 0.38, 0.38]);
skewmap(afrmask) = nan;
pcolor(ax4, longrid,latgrid,skewmap);
shading flat; 

cb4 = colorbar;
% set(cb4,'position',[0.05, 0.18, 0.02, 0.5])
title(cb4, 'Skewness', 'FontSize', 12)
caxis([-4 6]); 

% customize the colormap colors
GreenMap = [linspace(0.39, 1, 40)', linspace(0.83, 1, 40)', linspace(0.07, 1, 40)'];
PurpleMap = [linspace(1, 0.72, 60)', linspace(1, 0.27, 60)', linspace(1, 0.99, 60)'];
GWPMap = [GreenMap; PurpleMap];
colormap(ax4, GWPMap);

hold on

% Plot the African boundary
plot(coastline, 'EdgeColor', 'k', 'LineWidth', 2.5, 'FaceAlpha', 0);
% Plot Madgascar's Big Island
plot(madgline.Vertices(1:50:end, 1), madgline.Vertices(1:50:end, 2),'k', 'LineWidth', 2.5);
% Plot tectonic plate outlines
for ip = 1:117
    plot(Tplates(ip).X, Tplates(ip).Y, 'LineWidth', 2, 'Color',[0.4 0.4 0.4])
end

% Plot African cratons
plot(AC.Congo(:,2), AC.Congo(:,1), 'color', [0.2 0.2 0.2], 'lineWidth', 1);
plot(AC.Kala(:,2), AC.Kala(:,1), 'color', [0.2 0.2 0.2], 'lineWidth', 1);
plot(AC.Sahara(:,2), AC.Sahara(:,1), 'color', [0.2 0.2 0.2], 'lineWidth', 1);
plot(AC.Tanza(:,2), AC.Tanza(:,1), 'color', [0.2 0.2 0.2], 'lineWidth', 1);
plot(AC.West(:,2), AC.West(:,1), 'color', [0.2 0.2 0.2], 'lineWidth', 1);
plot(AC.Zim(:,2), AC.Zim(:,1), 'color', [0.2 0.2 0.2], 'lineWidth', 1);

xlim([minlon maxlon]);
ylim([minlat maxlat]);
axesm('mercator','maplatlimit',[minlat maxlat],'maplonlimit',[minlon maxlon]);

hold off


% set the size of the plot
x0=10;
y0=10;
fwidth=1000;
fheight=500;
set(gcf,'position',[x0,y0,fwidth,fheight]);

% save the figure
figpath = '/Users/sxue3/Documents/BayMap_Figures/fig/fig4_L35P_mapUncer.pdf';
set(gcf,'Units','Inches');
pos = get(gcf,'Position');
set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
print(gcf,figpath,'-painters', '-dpdf','-r0');