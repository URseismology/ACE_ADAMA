%%%% Fig 1: African topo map with hotspots and red shapes 
%%%% No stations will be marked in this figure and labels are added afterward
%%%% This code uses m_map
%%%% Siyu Xue -- Jan 10. 2023
addpath('/Users/sxue3/Documents/BayMap_Figures/Utils')
addpath('/Users/sxue3/Documents/ADAMA_Figures/m_map')
addpath('/Users/sxue3/Documents/Figures/m_map');

%% Read in data

% the hotspots coordinates
fileID = fopen('./Data/GeoData/hotspots.txt', 'r');
formatSpec = '%f %f';
A = textscan(fileID, formatSpec);
fclose(fileID);

hotlon = A{1,1};
hotlat = A{1,2};

load('/Users/sxue3/Documents/BayMap_Figures/Data/GeoData/AfricanCratons.mat'); % load cratons
load('/Users/sxue3/Documents/BayMap_Figures/Data/GeoData/AfricanCrustal.mat') % load crustal types

%% Plot the topo of Africa

fs = 18;  % font size
clc
clf

% set the project as a regular map plot
m_proj('lambert','lon',[-26 60],'lat',[-39 40]);
caxis([-6000 3000]);   % set the colorbar range
colormap([m_colmap('blues',200); m_colmap('gland', 100)])   % set the colorbar color

[ELEV,LONG,LAT] = m_etopo2([-26 60 -39 40]);    % get the elevation (topography)
m_image(LONG(1,:),LAT(:,1),ELEV);

% set the lat&lon grid and labels
m_grid('linestyle','none','tickdir','out','linewidth',3);

% set the colorbar display
topoax = m_contfbar(0.9,[.17 .79], ELEV',[-6000:3000],...
       'axfrac',.02,'endpiece','no','levels','match','edgecolor','none'); 
set(topoax,'fontsize',(fs-1));
ylabel(topoax, 'Topography (km)');

hold on;

% Plot the boundary lines and cratons
% -- boundary
% fname = '/Users/sxue3/Documents/Figures/geoData/AfricaCoast.mat';
% afcoast = load(fname);
% aflon = wrapTo180(afcoast.XY(:,1));
% aflat = afcoast.XY(:,2);
% m_line(aflon,aflat,'color','k','linewi',2);

% plot cratons
m_line(AC.Congo(:,2), AC.Congo(:,1), 'color',[0.3, 0.3, 0.3],'linewi',2);
m_line(AC.Kala(:,2), AC.Kala(:,1), 'color',[0.3, 0.3, 0.3],'linewi',2);
m_line(AC.Sahara(:,2), AC.Sahara(:,1), 'color',[0.3, 0.3, 0.3],'linewi',2);
m_line(AC.Tanza(:,2), AC.Tanza(:,1), 'color',[0.3, 0.3, 0.3],'linewi',2);
m_line(AC.West(:,2), AC.West(:,1), 'color',[0.3, 0.3, 0.3],'linewi',2);
m_line(AC.Zim(:,2), AC.Zim(:,1), 'color',[0.3, 0.3, 0.3],'linewi',2);

% -- plates
plotplates()


% -- Plot the hotspots
for i = 1:length(hotlon)
    lat = hotlat(i);
    lon = hotlon(i);
    m_line(lon, lat, 'marker','o', ...
        'color','w','linewi',1,...
        'linest','none','markersize',6,'markerfacecolor',[0.89, 0, 0.13]);
end

xs = 0:0.1:50;
ys1 = 24/25 * xs - 21.5;
m_line(xs, ys1, 'color', [0.502 0.502 0.502], 'linewi',2)


% plot crustal regions
m_line(GEO.MobileBelts.Other.Cap_belt(:,2), GEO.MobileBelts.Other.Cap_belt(:,1), 'color', 'k','linewi', 1)
m_line(GEO.MobileBelts.Other.Kibaran_belt(:,2), GEO.MobileBelts.Other.Kibaran_belt(:,1), 'color', 'k','linewi', 1)
m_line(GEO.MobileBelts.Other.Namaquanatal(:,2), GEO.MobileBelts.Other.Namaquanatal(:,1), 'color', 'k','linewi', 1)
m_line(GEO.MobileBelts.Other.Ruwenzky(:,2), GEO.MobileBelts.Other.Ruwenzky(:,1), 'color', 'k','linewi', 1)
m_line(GEO.MobileBelts.Other.Mauritania(:,2), GEO.MobileBelts.Other.Mauritania(:,1), 'color', 'k','linewi', 1)
m_line(GEO.MobileBelts.WestAfrica_belt(:,2), GEO.MobileBelts.WestAfrica_belt(:,1), 'color', 'k','linewi', 1)
m_line(GEO.MobileBelts.Ouban_Damara.Damara(:,2), GEO.MobileBelts.Ouban_Damara.Damara(:,1), 'color', 'k','linewi', 1)
m_line(GEO.MobileBelts.Ouban_Damara.Oubangides(:,2), GEO.MobileBelts.Ouban_Damara.Oubangides(:,1), 'color', 'k','linewi', 1)

m_line(GEO.Orogens.AtlasMountains(:,2), GEO.Orogens.AtlasMountains(:,1), 'color', 'k','linewi', 1)
m_line(GEO.Orogens.EastAfrica.EastAfrica(:,2), GEO.Orogens.EastAfrica.EastAfrica(:,1), 'color', 'k','linewi', 1)
m_line(GEO.Orogens.EastAfrica.NorthMozamb(:,2), GEO.Orogens.EastAfrica.NorthMozamb(:,1), 'color', 'k','linewi', 1)

m_line(GEO.Basins.Congo(:,2), GEO.Basins.Congo(:,1), 'color', 'k','linewi', 1)
m_line(GEO.Basins.Tindouf.Tindouf(:,2), GEO.Basins.Tindouf.Tindouf(:,1), 'color', 'k','linewi', 1)
m_line(GEO.Basins.Tindouf.MauriSenba(:,2), GEO.Basins.Tindouf.MauriSenba(:,1), 'color', 'k','linewi', 1)
m_line(GEO.Basins.Taoudeni.B1(:,2), GEO.Basins.Taoudeni.B1(:,1), 'color', 'k','linewi', 1)

m_line(GEO.Archean.WestAfrica.LeoManShield(:,2), GEO.Archean.WestAfrica.LeoManShield(:,1), 'color', 'k','linewi', 1)
m_line(GEO.Archean.WestAfrica.Reguibatshield(:,2), GEO.Archean.WestAfrica.Reguibatshield(:,1), 'color', 'k','linewi', 1)
m_line(GEO.Archean.Tanzania(:,2), GEO.Archean.Tanzania(:,1), 'color', 'k','linewi', 1)
m_line(GEO.Archean.CongoAll.AC(:,2), GEO.Archean.CongoAll.AC(:,1), 'color', 'k','linewi', 1)
m_line(GEO.Archean.CongoAll.Bomu(:,2), GEO.Archean.CongoAll.Bomu(:,1), 'color', 'k','linewi', 1)
% plot(GEO.Archean.CongoAll.Congo(:,2), GEO.Archean.CongoAll.Congo(:,1), 'k','lineWidth', 1)
m_line(GEO.Archean.CongoAll.Kasai(:,2), GEO.Archean.CongoAll.Kasai(:,1), 'color', 'k','linewi', 1)
m_line(GEO.Archean.CongoAll.Ntem(:,2), GEO.Archean.CongoAll.Ntem(:,1), 'color', 'k','linewi', 1)
m_line(GEO.Archean.SaharaMeta.A1(:,2), GEO.Archean.SaharaMeta.A1(:,1), 'color', 'k','linewi', 1)
m_line(GEO.Archean.SaharaMeta.A2(:,2), GEO.Archean.SaharaMeta.A2(:,1), 'color', 'k','linewi', 1)
m_line(GEO.Archean.SaharaMeta.A3(:,2), GEO.Archean.SaharaMeta.A3(:,1), 'color', 'k','linewi', 1)
m_line(GEO.Archean.SaharaMeta.A4(:,2), GEO.Archean.SaharaMeta.A4(:,1), 'color', 'k','linewi', 1)
m_line(GEO.Archean.WestAfrica_MZ.Dahomeyshield(:,2), GEO.Archean.WestAfrica_MZ.Dahomeyshield(:,1), 'color', 'k','linewi', 1)
m_line(GEO.Archean.WestAfrica_MZ.Tuaregshield(:,2), GEO.Archean.WestAfrica_MZ.Tuaregshield(:,1), 'color', 'k','linewi', 1)
m_line(GEO.Archean.KaapVaal(:,2), GEO.Archean.KaapVaal(:,1), 'color', 'k','linewi', 1)




% set the size of the plot
x0=10;
y0=10;
width=1200;
height=800;
set(gcf,'position',[x0,y0,width,height]);

% save the figure
fig = gcf;
saveFig('fig1_AfrTopo.pdf', '/Users/sxue3/Documents/BayMap_Figures/fig/', 1, fig);
