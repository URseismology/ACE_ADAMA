%%%% Fig 1: African map with station marked on them
%%%% This code uses m_map
%%%% Input files: ppconnCnt.csv and African craton file
%%%% Siyu Xue -- Jan 10. 2023
addpath('/Users/sxue3/Documents/BayMap_Figures/Utils')

%% Read in Data
dat = readtable('/Users/sxue3/Documents/ADAMA_Figures/data/ppconnCnt.csv');
stas = dat.Station;
nets = dat.Network;
nSta = length(stas);
lats = dat.Latitude;
lons = dat.Longitude;
connCnt = dat.connCnt;

S1 = load('/Users/sxue3/Documents/BayMap_Figures/Data/AfrCountry.mat').S1;
%% Plot the boundary of Africa and stations
addpath('/Users/sxue3/Documents/ADAMA_figures/m_map');
fs = 18;  % font size
clc

% plot map with topo
clf
hold on;

m_proj('lambert','lon',[-26 60],'lat',[-39 40]);  % African boundaries
m_grid('linestyle','none','tickdir','out','linewidth',3);


% Plot the country lines and cratons
% -- countries
for iCntry = 1:length(S1)
    CX = [S1(iCntry).X];
    CY = [S1(iCntry).Y];

    m_line(CX(1:10:end), CY(1:10:end), 'color', [0.5 0.5 .5], 'linewi',2);  % states ...
end

% -- cratons
load('/Users/sxue3/Documents/BayMap_Figures/Data/AfricanCratons.mat');
% plot cratons
m_line(AC.Congo(:,2), AC.Congo(:,1), 'color',[0.3, 0.3, 0.3],'linewi',2);
m_line(AC.Kala(:,2), AC.Kala(:,1), 'color',[0.3, 0.3, 0.3],'linewi',2);
m_line(AC.Sahara(:,2), AC.Sahara(:,1), 'color',[0.3, 0.3, 0.3],'linewi',2);
m_line(AC.Tanza(:,2), AC.Tanza(:,1), 'color',[0.3, 0.3, 0.3],'linewi',2);
m_line(AC.West(:,2), AC.West(:,1), 'color',[0.3, 0.3, 0.3],'linewi',2);
m_line(AC.Zim(:,2), AC.Zim(:,1), 'color',[0.3, 0.3, 0.3],'linewi',2);

% % -- plot grey lines (cross-sections)
xs = 0:0.1:50;
ys1 = 24/25 * xs - 21.5;
m_line(xs, ys1, 'color', [0.57 0.36 0.51], 'linewi',2)

% -- set the color of each station marker
for iSta = 1 : nSta
    lat = lats(iSta);
    lon = lons(iSta);
        
    m_line(lon, lat, 'marker','o', 'color','r','linewi',1,...
        'linest','none','markersize',4);
end


% set the size of the plot
x0=10;
y0=10;
width=800;
height=800;
set(gcf,'position',[x0,y0,width,height]);

% save the figure
fig = gcf;
% saveFig('fig1_AfrSta.pdf', '/Users/sxue3/Documents/BayMap_Figures/fig/', 1, fig);
