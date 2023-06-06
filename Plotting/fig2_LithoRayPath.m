%%% Plot the Litho RayPath coverage in Africa (not included in the paper)
%%% Top station are marked using special markers
%%% This code uses m_map
%%% Siyu Xue -- Feb 10. 2023
addpath('/Users/sxue3/Documents/BayMap_Figures/Utils')
addpath('/Users/sxue3/Documents/ADAMA_Figures/m_map');

%% Read in Data
fileID = fopen('./Data/LithoRP.txt', 'r');
formatSpec = '%f %f %f %f';
A = textscan(fileID, formatSpec);
fclose(fileID);
 
lithlat = A{1,1};
lithlon = A{1,2};
lithcount = A{1,4};

logcount = log10(lithcount);

S1 = load('/Users/sxue3/Documents/BayMap_Figures/Data/AfrCountry.mat').S1;
madgline = polyshape((S1(448).X).',(S1(448).Y).');  % load Madgascar's Big Island

fname = '/Users/sxue3/Documents/BayMap_Figures/Data/AfricaCoast.mat';
afcoast = load(fname);
aflon = wrapTo180(afcoast.XY(:,1));
aflat = afcoast.XY(:,2);
coastline = polyshape(aflon, aflat);

fs = 18;  % font size
% -- cratons
load('/Users/sxue3/Documents/ADAMA_Figures/data/cratons/AfricanCratons.mat');

% interpolate the Litho coverage map
F = scatteredInterpolant(lithlat,lithlon,logcount, 'natural');
[LG,LT]=meshgrid(-26:60,-39:40);  % interpolate to one degree
InterpDat = F(LG, LT);

% Create a mask to seperate region within Africa
afrmask = zeros(80, 87);
mapsize = size(afrmask);
for xi = 1:mapsize(2)
    loni = LG(1,xi);
    for yi = 1:mapsize(1)
        lati = LT(yi,1);

        if isinterior(coastline, loni, lati) == 1 || isinterior(madgline, loni, lati) == 1 ...
            || isinterior(coastline, loni+1, lati+1) == 1 || isinterior(madgline, loni+1, lati+1) == 1
            afrmask(yi, xi) = 1;
        end
    end
end

afrmask = logical(afrmask - 1);  % points not in Africa == 1
InterpDat(afrmask) = nan;

%% Plot the Litho RayPath of Africa
hold on
m_proj('lambert','lon',[-26 60],'lat',[-39 40]);  % African boundaries
m_pcolor(LG,LT,InterpDat)
m_grid('linestyle','none','tickdir','out','linewidth',3);

colormap(flipud(jet));
colorbar('fontsize',(fs-1));
caxis([0.5 4.1]);

% -- boundary
m_line(aflon,aflat,'color','k','linewi',2);

% plot Madagascar
Malon = (S1(448).X).';
Malat = (S1(448).Y).';
Malon = Malon(1:10:end);
Malat = Malat(1:10:end);
m_line(Malon ,Malat,'color','k','linewi',2);


% plot cratons
m_line(AC.Congo(:,2), AC.Congo(:,1), 'color','k','linewi',2);
m_line(AC.Kala(:,2), AC.Kala(:,1), 'color','k','linewi',2);
m_line(AC.Sahara(:,2), AC.Sahara(:,1), 'color','k','linewi',2);
m_line(AC.Tanza(:,2), AC.Tanza(:,1), 'color','k','linewi',2);
m_line(AC.West(:,2), AC.West(:,1), 'color','k','linewi',2);
m_line(AC.Zim(:,2), AC.Zim(:,1), 'color','k','linewi',2);


% set the size of the plot
x0=10;
y0=10;
width=1500;
height=700;
set(gcf,'position',[x0,y0,width,height]);

% %% Plot the difference between ADAMA and Litho
% ADAMAF = load('/Users/sxue3/Documents/BayMap_Figures/Data/ADAMA_raypath.mat').F;
% ADval = ADAMAF(LG, LT);
% ADval(isnan(ADval))=0;
% ADval(afrmask) = nan;
% 
% Diffvalue = ADval./(10.^InterpDat);
% 
% figure(24)
% m_proj('lambert','lon',[-26 60],'lat',[-39 40]);  % African boundaries
% m_pcolor(LG,LT,Diffvalue)
% m_grid('linestyle','none','tickdir','out','linewidth',3);
% % -- boundary
% m_line(aflon,aflat,'color','k','linewi',2);
% 
% colormap(m_colmap('jet'));
% colorbar('fontsize',(fs-1));
% title('ADAMA Times of Litho')
% 
% % set the size of the plot
% x0=10;
% y0=10;
% width=1500;
% height=700;
% set(gcf,'position',[x0,y0,width,height]);


% save the figure
% fig = gcf;
% saveFig('fig2_AfrStaQlty_pp.pdf', '/Users/sxue3/Documents/Figures/fig/', 1, fig);