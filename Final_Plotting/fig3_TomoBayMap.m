%%%% Fig 3: plot a Bayesian tomo-map snapshot and other statistics
%%%%% Siyu Xue -- Feb 10. 2023
addpath('/Users/sxue3/Documents/BayMap_Figures/Utils')
addpath('/Users/sxue3/Documents/ADAMA_Figures')

%% Read in data
part_xctrs = load('/Users/sxue3/Documents/BayMap_Figures/Data/LovePData/L35/xctrs.mat').xctrs;
part_yctrs = load('/Users/sxue3/Documents/BayMap_Figures/Data/LovePData/L35/yctrs.mat').yctrs;
part_vels = load('/Users/sxue3/Documents/BayMap_Figures/Data/LovePData/L35/velvals.mat').velvals;
longrid = load('/Users/sxue3/Documents/BayMap_Figures/Data/RayleighPData/plotpara.mat').longrid;
latgrid = load('/Users/sxue3/Documents/BayMap_Figures/Data/RayleighPData/plotpara.mat').latgrid;

% load African country boundaries
S1 = load('/Users/sxue3/Documents/BayMap_Figures/Data/AfrCountry.mat').S1;
madgline = polyshape((S1(448).X).',(S1(448).Y).');  % load Madgascar's Big Island

% load the coast
afcoast = load('/Users/sxue3/Documents/BayMap_Figures/Data/AfricaCoast.mat');
aflon = wrapTo180(afcoast.XY(:, 1));
aflat = afcoast.XY(:,2);
coastline = polyshape(aflon, aflat);

% load cratons
load('/Users/sxue3/Documents/BayMap_Figures/Data/AfricanCratons.mat');

% load the Tectonic Plate Boundaries
Tplates = shaperead('./Data/Tectonic_PB.shp');

% Compute the velocity final map
nummaps = size(part_vels);
i = nummaps(1);

% x coordinates
x_thisiter = part_xctrs(i,:);
x_thisiter=x_thisiter.';
% y coordinates
y_thisiter = part_yctrs(i,:);
y_thisiter=y_thisiter.';
% Velocities
vel_thisiter = part_vels(i,:);
vel_thisiter=vel_thisiter.';

% Interpolate!
velmap_thisiter = interp_voronoi(x_thisiter,y_thisiter,vel_thisiter,...
longrid,latgrid);

[vx, vy] = voronoi(x_thisiter, y_thisiter);

minlat = min(latgrid(:,1)); maxlat = max(latgrid(:,1));
minlon = min(longrid(1,:)); maxlon = max(longrid(1,:));


% If plotting using CellVelDistribution data (then can run the 'Plots Aggregation' section directly)
cellcenters = zeros(60,65);
x_vector = linspace(minlon, maxlon, 65);
y_vector = linspace(minlat, maxlat, 60);
[Xlon,Ylat] = meshgrid(x_vector,y_vector);
nummaps = nummaps(1);
load('Data/LCellVelDist/L35_P_dist.mat')

%% Compute the density map

% Dimensions: Longitude x latitude 
cellcenters = zeros(60,65);
x_vector = linspace(minlon, maxlon, 65);
y_vector = linspace(minlat, maxlat, 60);

% Iterate through each velocity map, and add its cell center to the total
nummaps = size(part_vels);
nummaps = nummaps(1);

totsum = 0;
for i = nummaps/2: nummaps%1:nummaps
    totsum = totsum + 1;
    %disp([totsum, i]);
    
    % x coordinates
    x_thisiter = part_xctrs(i,:);
    x_thisiter = x_thisiter.';
    % y coordinates
    y_thisiter = part_yctrs(i,:);
    y_thisiter = y_thisiter.';
    % Velocities
    vel_thisiter = part_vels(i,:);
    vel_thisiter = vel_thisiter.';


    for icell = 1:length(x_thisiter)

        if x_thisiter(icell) == 0 && y_thisiter(icell) == 0
            continue
        end

        if x_thisiter(icell) < minlon || x_thisiter(icell) > maxlon
            continue
        end

        if y_thisiter(icell) < minlat || y_thisiter(icell) > maxlat
            continue
        end

        x_grid = find(x_vector <= x_thisiter(icell), 1, 'last');
        y_grid = find(y_vector <= y_thisiter(icell), 1, 'last');

        if vel_thisiter(icell) > 5  % interested in the density of cells where vel > 4.5
        cellcenters(y_grid,x_grid) = cellcenters(y_grid,x_grid) + 1; % add a cell count to the grid
        end
    end    
end

densmap = cellcenters ./ totsum;  % compute the density
[Xlon,Ylat] = meshgrid(x_vector,y_vector);


%% Compute total points in Africa over iterations
% Iterate through each velocity map, and add its cell center to the total
nummaps = size(part_vels);

cellinafr = zeros(nummaps/100, 1);
isample = 0;
for i = 1:100:nummaps

    % x coordinates
    x_thisiter = part_xctrs(i,:);
    x_thisiter = x_thisiter.';
    % y coordinates
    y_thisiter = part_yctrs(i,:);
    y_thisiter = y_thisiter.';

    iruncell = 0;
    for icell = 1:length(x_thisiter)
        x = x_thisiter(icell);
        y = y_thisiter(icell);
        if x == 0 && y == 0
            continue
        end

        if isinterior(coastline, x, y) == 1 || isinterior(madgline, x, y) == 1
            iruncell = iruncell + 1;
        end
    end
    isample = isample + 1;
    cellinafr(isample) = iruncell;
end

%% Collect velocity in Africa all iterations

% Iterate through each velocity map, and add its cell center to the total
nummaps = size(part_vels, 1);
nummaps = nummaps(1);

velinafr = [];

for i = 1:100:nummaps
    
    % x coordinates
    x_thisiter = part_xctrs(i,:);
    x_thisiter = x_thisiter.';
    % y coordinates
    y_thisiter = part_yctrs(i,:);
    y_thisiter = y_thisiter.';
    % Velocities
    vel_thisiter = part_vels(i,:);
    vel_thisiter=vel_thisiter.';

    for icell = 1:length(x_thisiter)

        x = x_thisiter(icell);
        y = y_thisiter(icell);
        if x == 0 && y == 0
            continue
        end

        if isinterior(coastline, x, y) == 1 || isinterior(madgline, x, y) == 1
            velinafr = [velinafr, vel_thisiter(icell)];
        end
    end
end

%% Plots Aggregation 
set(0,'defaultfigurecolor',[1 1 1])

% (a1) Plot the velocity final map
ax1 = axes();
avgplot=pcolor(ax1,longrid,latgrid,velmap_thisiter);  % plot cells and values
shading flat; 
% caxis([1 5]); cb1 = colorbar;
caxis([3.8 5]); cb1 = colorbar;
annotation('textbox',[.83 .58 .3 .3], 'String','Velocity (km/s)','EdgeColor','none', 'FontSize', 12)
% title(cb1, 'Velocity (m/s)', 'FontSize', 12)

% customize the colormap colors
RedMap = [linspace(0.98, 1, 128)', linspace(0.16, 1, 128)', linspace(0.16, 1, 128)'];
BlueMap = [linspace(1, 0.16, 128)', linspace(1, 0.57, 128)', linspace(1, 0.98, 128)'];
RWBMap = [RedMap; BlueMap];
colormap(ax1, RWBMap);

hold on
plot(vx, vy, 'k');  % plot the cell boundary

% plot the cell centers (if in Africa)
inAfr = ones(length(x_thisiter), 1);  % indicate if a point is in Africa & Madg
for icell = 1:length(x_thisiter)
    x = x_thisiter(icell);
    y = y_thisiter(icell);
    if isinterior(coastline, x, y) == 0 && isinterior(madgline, x, y) == 0
        inAfr(icell) = 0;
    end
end

plot(x_thisiter(logical(inAfr)), y_thisiter(logical(inAfr)), 'k.', 'MarkerSize', 6)  % plot the center of each cell
plot(coastline, 'EdgeColor', 'w', 'LineWidth', 5, 'FaceAlpha', 0);  % plot the African boundary

% Plot Madgascar's Big Island
plot(madgline.Vertices(1:50:end, 1), madgline.Vertices(1:50:end, 2),'w', 'LineWidth', 5);

% Plot tectonic plate outlines
for ip = 1:117
    plot(Tplates(ip).X, Tplates(ip).Y, 'k', 'LineWidth', 3)
end

% Plot African cratons
plot(AC.Congo(:,2), AC.Congo(:,1), 'color', [0.4 0.4 0.4], 'lineWidth', 2.5);
plot(AC.Kala(:,2), AC.Kala(:,1), 'color', [0.4 0.4 0.4], 'lineWidth', 2.5);
plot(AC.Sahara(:,2), AC.Sahara(:,1), 'color', [0.4 0.4 0.4], 'lineWidth', 2.5);
plot(AC.Tanza(:,2), AC.Tanza(:,1), 'color', [0.4 0.4 0.4], 'lineWidth', 2.5);
plot(AC.West(:,2), AC.West(:,1), 'color', [0.4 0.4 0.4], 'lineWidth', 2.5);
plot(AC.Zim(:,2), AC.Zim(:,1), 'color', [0.4 0.4 0.4], 'lineWidth', 2.5);

xlim([minlon maxlon]);
ylim([minlat maxlat]);
axesm('mercator','maplatlimit',[minlat maxlat],'maplonlimit',...
  [minlon maxlon]);
xticks(10*ceil(minlon/10):10:10*floor(maxlon/10));
yticks(10*ceil(minlat/10):10:10*floor(maxlat/10));

annotation('textbox',[.125 .82 .075 .03], 'String', ...
    'L35-P','EdgeColor','k', 'BackgroundColor','w', 'FontSize', 18)
hold off

% (a2) Plot velocity in Africa all iterations
% [xstart ystart xend-xstart yend-ystart]
ax2bg = axes('position',[0.658, 0.73, 0.2, 0.18]);  % add white background for (a2)
ax2bg.YAxis.Visible = 'off'; % remove y-axis
ax2bg.XAxis.Visible = 'off'; % remove x-axis

ax2 = axes('position',[0.685, 0.75, 0.17, 0.10]);
hold on
h100 = histogram(velinafr.', 50);
% xlabel('Velocity (m/s)')
yticks(linspace(0, max(h100.Values), 6))
yticklabels({'0','0.2','0.4','0.6','0.8','1'})
xlim([0, 6])
alpha(.5)

totalvel = length(velinafr);
h50 = histogram(velinafr(1:0.5*totalvel).');
alpha(.5)
h10 = histogram(velinafr(1:0.1*totalvel).');
alpha(.5)

legend('all iter', 'first 50% iter', 'first 10% iter', 'Location',[0.56, 0.812, 0.08, 0.04]);


% (b1) Plot density map
% [xstart ystart xend-xstart yend-ystart]
ax3bg = axes('position',[0.1234, 0.1823, 0.28, 0.34]);  % add white background for (b1&b2)
% ax3bg = axes('position',[0.1234, 0.1823, 0.28, 0.346]);  % add white background for (b1&b2)
ax3bg.YAxis.Visible = 'off'; % remove y-axis
ax3bg.XAxis.Visible = 'off'; % remove x-axis

ax3 = axes('position',[0.125, 0.155, 0.27, 0.3]);
denplot = pcolor(ax3, Xlon,Ylat,densmap);
cb2 = colorbar();
title(cb2, 'Density', 'FontSize', 12)
colormap(ax3, hot);
set(cb2,'position',[0.92, 0.18, 0.02, 0.5])

hold on
% Plot African countries
for iCntry = 1:length(S1)
    CX = [S1(iCntry).X];
    CY = [S1(iCntry).Y];

    plot(ax3, CX(1:10:end), CY(1:10:end), 'color', [0.5 0.5 .5], 'linewi',1);  % states ...
end

% Plot the African boundary
plot(ax3, coastline, 'EdgeColor', 'w', 'LineWidth', 3, 'FaceAlpha', 0);  
% Plot Madgascar's Big Island
plot(ax3, madgline.Vertices(1:50:end, 1), madgline.Vertices(1:50:end, 2),'w', 'LineWidth', 3);

% Plot tectonic plate outlines
for ip = 1:117
    plot(Tplates(ip).X, Tplates(ip).Y, 'Color', '[0.6 0.6 0.6]', 'LineWidth', 2)
end

% Plot African cratons
plot(ax3,AC.Congo(:,2), AC.Congo(:,1), 'color', [0.85 0.85 0.85], 'lineWidth', 1.5);
plot(ax3,AC.Kala(:,2), AC.Kala(:,1), 'color', [0.85 0.85 0.85], 'lineWidth', 1.5);
plot(ax3,AC.Sahara(:,2), AC.Sahara(:,1), 'color', [0.85 0.85 0.85], 'lineWidth', 1.5);
plot(ax3,AC.Tanza(:,2), AC.Tanza(:,1), 'color', [0.85 0.85 0.85], 'lineWidth', 1.5);
plot(ax3,AC.West(:,2), AC.West(:,1), 'color', [0.85 0.85 0.85], 'lineWidth', 1.5);
plot(ax3,AC.Zim(:,2), AC.Zim(:,1), 'color', [0.85 0.85 0.85], 'lineWidth', 1.5);

xlim([minlon maxlon]);
ylim([minlat maxlat]);
axesm('mercator','maplatlimit',[minlat maxlat],'maplonlimit',...
  [minlon maxlon]);
hold off

% (b2) Plot total points in Africa over iterations
ax4 = axes('position',[0.155, 0.45, 0.24, 0.07]);

plot(cellinafr, 'Color', [0.6 0.6 0.6]);
hold on
plot(smooth(cellinafr, 20), 'b', 'LineWidth', 1.5)
xticks(linspace(0, nummaps/100, 6))
xticklabels({'0','0.2','0.4','0.6','0.8','1'})
xlim([0, nummaps/100])
annotation('textbox',[.155 .22 .3 .3], 'String', ...
    [num2str(nummaps), ' iterations'],'EdgeColor','none', 'FontSize', 12)
annotation('textbox',[.123 .522 .12 .02], 'String', ...
    ['cell count: ', num2str(sum(x_thisiter(logical(inAfr))~=0))],'EdgeColor','none', 'BackgroundColor','w', 'FontSize', 12)
% annotation('textbox',[.119 .526 .12 .018], 'String', ...
%     ['cell count: ', num2str(length(x_thisiter))],'EdgeColor','none', 'BackgroundColor','w', 'FontSize', 12)

% set the size of the plot
x0=10;
y0=10;
fwidth=800;
fheight=800;
set(gcf,'position',[x0,y0,fwidth,fheight]);

% save the figure
figpath = '/Users/sxue3/Documents/BayMap_Figures/fig/fig3_sum_L35P.pdf';
set(gcf,'Units','Inches');
pos = get(gcf,'Position');
set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
print(gcf,figpath,'-painters', '-dpdf','-r0');

% save('./R40_G.mat', 'densmap', 'cellinafr', 'velinafr');

