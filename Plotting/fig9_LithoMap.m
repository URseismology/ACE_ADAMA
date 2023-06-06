% Twicked from Yuri's plotting code 
% Plot Average Velocity map from Litho1.0 Map
% Input data is generated from Prj5/4_Bin/Yuri_plotAfrica/plotLITHOmap_Ma2014.m
% Nov. 30, 2021

%% Load data
addpath('/Users/sxue3/Documents/BayMap_Figures/Utils')

% load African country boundaries
S1 = load('/Users/sxue3/Documents/BayMap_Figures/Data/AfrCountry.mat').S1;
madgline = polyshape((S1(448).X).',(S1(448).Y).');  % load Madgascar's Big Island

% load the coast
afcoast = load('/Users/sxue3/Documents/BayMap_Figures/Data/AfricaCoast.mat');
aflon = wrapTo180(afcoast.XY(:, 1));
aflat = afcoast.XY(:,2);
coastline = polyshape(aflon, aflat);

% load cratons
AC = load('/Users/sxue3/Documents/ADAMA_Figures/data/cratons/AfricanCratons.mat');
latCongo = AC.Congo(:,1);
lonCongo = AC.Congo(:,2);
latKala = AC.Kala(:,1);
lonKala = AC.Kala(:,2);
latSahara = AC.Sahara(:,1);
lonSahara = AC.Sahara(:,2);
latTanza = AC.Tanza(:,1);
lonTanza = AC.Tanza(:,2);
latWest = AC.West(:,1);
lonWest = AC.West(:,2);
latZim = AC.Zim(:,1);
lonZim = AC.Zim(:,2);

%% Plot phase velocity map
lithofig=figure();
lithofig.Units='normalized';

lithofig.OuterPosition(3)=0.4;
lithofig.OuterPosition(4)=0.85;

minlon = min(longrid(1,:));
minlat = min(latgrid(:,1));
maxlon = max(longrid(1,:));
maxlat = max(latgrid(:,1));

axesm('mercator','maplatlimit',[minlat maxlat],'maplonlimit',...
  [minlon maxlon])

% Axes and title
ax=gca;
ax.XLim=[minlon maxlon];
ax.YLim=[minlat maxlat];

ax.XLabel.String='Longitude (Degrees)';
ax.XLabel.FontSize=13;
ax.YLabel.String='Latitude (Degrees)';
ax.YLabel.FontSize=13;

% Convert frequency in mHz to period
pdval=period;
% titlestr1=sprintf('%s Wave Phase Velocities (T = %.2fs)',wavestr,pdval);
% titlestr2='Ma et al. (2014) Dispersion Model (1 x 1 Degree)';
% ax.Title.String={titlestr1;titlestr2};

ax.XTick=10*ceil(minlon/10):10:10*floor(maxlon/10);
ax.YTick=10*ceil(minlat/10):10:10*floor(maxlat/10);

ax.FontSize=12;
ax.Title.FontSize=12;
hold on

% customize the colormap colors
RedMap = [linspace(0.98, 1, 128)', linspace(0.16, 1, 128)', linspace(0.16, 1, 128)'];
BlueMap = [linspace(1, 0.16, 128)', linspace(1, 0.57, 128)', linspace(1, 0.98, 128)'];
RWBMap = [RedMap; BlueMap];
colormap(ax, RWBMap)
cbar=colorbar;
cbar.Label.String='Phase Velocity (km/s)';
caxis([3.55 4.7]);   % set the range of color bar manuallty
cbar.FontSize=12;

% Plot
lithoplot=pcolor(longrid,latgrid,velgrid);
lithoplot.EdgeAlpha=0;

% Plot African boundary
plot(coastline, 'EdgeColor', 'k', 'LineWidth', 5, 'FaceAlpha', 0);  % plot the African boundary
% Plot Madgascar's Big Island
plot(madgline.Vertices(1:50:end, 1), madgline.Vertices(1:50:end, 2),'k', 'LineWidth', 5);
% Plot tectonic plate outlines
plotplates();

plot(lonCongo, latCongo, 'color', 'k', 'lineWidth', 2);
plot(lonKala, latKala, 'color', 'k', 'lineWidth', 2);
plot(lonSahara, latSahara, 'color', 'k', 'lineWidth', 2);
plot(lonTanza, latTanza, 'color', 'k', 'lineWidth', 2);
plot(lonWest, latWest, 'color', 'k', 'lineWidth', 2);
plot(lonZim, latZim, 'color', 'k', 'lineWidth', 2);

% Axes limits
ax.XLim=[-20 55];
ax.YLim=[-38 40];

%% save the figure
% saveFig('fig9_LithoMap.pdf', '/Users/sxue3/Documents/Figures/fig/', 1, gcf);