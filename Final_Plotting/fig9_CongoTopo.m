%%%%% Fig 9: plot the African topo map, zoom into the transec
%%%%% No stations will be marked in this figure
%%%%% This code uses m_map
%%%%% Siyu Xue -- May 10. 2023

addpath('/Users/sxue3/Documents/BayMap_Figures/Utils')
addpath('/Users/sxue3/Documents/ADAMA_Figures/m_map');
clear

%% Plot the topo of Africa
clear

fs = 18;  % font size
clc

% plot map with topo
clf
% set the project as oblique (projection)
m_proj('oblique mercator','longitudes',[43 10], ...
           'latitudes',[20 -12],'direction','vertical','aspect',.5);

caxis([-6000 3000]);    % set the colorbar range
colormap([m_colmap('blues',200); m_colmap('gland', 100)])   % set the colorbar color

[ELEV,LONG,LAT] = m_etopo2([0 55 -25 30]);  % get the elevation (topography)
m_image(LONG(1,:),LAT(:,1),ELEV);

% set the lat&lon grid and labels
% m_grid('xlabeldir','end','xtick',8,'ytick',8,'fontsize',10,...
%     'xaxisLocation','bottom', 'yaxisLocation','right');
m_grid('tickdir','out');

% set the colorbar display
topoax = m_contfbar([.1 .9], 0 ...
    , ELEV',[-3000:3000],...
       'axfrac',.02,'endpiece','no','levels','match','edgecolor','none'); 
set(topoax,'fontsize',(fs-1));
title(topoax, 'Topography (km)');

hold on;

% -- cratons
load('/Users/sxue3/Documents/BayMap_Figures/Data/AfricanCratons.mat');
% plot cratons
m_line(AC.Congo(:,2), AC.Congo(:,1), 'color','k','linewi',2);
m_line(AC.Kala(:,2), AC.Kala(:,1), 'color','k','linewi',2);
m_line(AC.Sahara(:,2), AC.Sahara(:,1), 'color','k','linewi',2);
m_line(AC.Tanza(:,2), AC.Tanza(:,1), 'color','k','linewi',2);
m_line(AC.West(:,2), AC.West(:,1), 'color','k','linewi',2);
m_line(AC.Zim(:,2), AC.Zim(:,1), 'color','k','linewi',2);

m_line(AC.Ntem_craton(:,2), AC.Ntem_craton(:,1), 'color',[0.3, 0.3, 0.3],'linewi',2);
m_line(AC.Bomu_craton(:,2), AC.Bomu_craton(:,1), 'color',[0.3, 0.3, 0.3],'linewi',2);
m_line(AC.Kasai_craton(:,2), AC.Kasai_craton(:,1), 'color',[0.3, 0.3, 0.3],'linewi',2);
m_line(AC.Angola_craton(:,2), AC.Angola_craton(:,1), 'color',[0.3, 0.3, 0.3],'linewi',2);
m_line(AC.Congo_basin(:,2), AC.Congo_basin(:,1), 'color',[0.57 0.36 0.51],'linewi',2);

% -- plates
plotplates()

xs = 10:0.1:45;
ys1 = 24/25 * xs - 21.5;
m_line(xs, ys1, 'color', [1 0 0], 'linewi',2)

% add the north arrow
m_northarrow(35,22,4,'type',2);

% set the size of the plot
x0=50;
y0=50;
width=300;
height=700;
set(gcf,'position',[x0,y0,width,height]);

% save the figure
fig = gcf;
% saveFig('fig9_Congotopo.pdf', '/Users/sxue3/Documents/BayMap_Figures/fig/', 1, fig);
