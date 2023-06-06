%%%% Fig 1: African topo map with hotspots and red shapes 
%%%% No stations will be marked in this figure and labels are added afterward
%%%% This code uses m_map
%%%% Siyu Xue -- Jan 10. 2023
addpath('/Users/sxue3/Documents/BayMap_Figures/Utils')

%% Read in data

% the hotspots coordinates
fileID = fopen('./Data/hotspots.txt', 'r');
formatSpec = '%f %f';
A = textscan(fileID, formatSpec);
fclose(fileID);

hotlon = A{1,1};
hotlat = A{1,2};

%% Plot the topo of Africa
addpath('/Users/sxue3/Documents/Figures/m_map');
fs = 18;  % font size
clc

% plot map with topo
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

% -- cratons
load('/Users/sxue3/Documents/BayMap_Figures/Data/AfricanCratons.mat');
% plot cratons
m_line(AC.Congo(:,2), AC.Congo(:,1), 'color',[0.3, 0.3, 0.3],'linewi',2);
m_line(AC.Kala(:,2), AC.Kala(:,1), 'color',[0.3, 0.3, 0.3],'linewi',2);
m_line(AC.Sahara(:,2), AC.Sahara(:,1), 'color',[0.3, 0.3, 0.3],'linewi',2);
m_line(AC.Tanza(:,2), AC.Tanza(:,1), 'color',[0.3, 0.3, 0.3],'linewi',2);
m_line(AC.West(:,2), AC.West(:,1), 'color',[0.3, 0.3, 0.3],'linewi',2);
m_line(AC.Zim(:,2), AC.Zim(:,1), 'color',[0.3, 0.3, 0.3],'linewi',2);

% -- plates
plotplates()


% % -- plot all the red "cycle" 
% for i=1:60
%     fileID = fopen(['./Data/red_shapes/', num2str(i),'.txt'], 'r');
%     formatSpec = '%f %f';
%     A = textscan(fileID, formatSpec);
%     fclose(fileID);
% 
%     redlon = A{1,1};
%     redlat = A{1,2};
% 
%     % plot the red areas
%     m_line(redlon, redlat, 'color', [0.89, 0, 0.13], 'linewi', 3);
% end

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

% set the size of the plot
x0=10;
y0=10;
width=1200;
height=800;
set(gcf,'position',[x0,y0,width,height]);

% save the figure
fig = gcf;
% saveFig('fig1_AfrTopo.pdf', '/Users/sxue3/Documents/BayMap_Figures/fig/', 1, fig);
