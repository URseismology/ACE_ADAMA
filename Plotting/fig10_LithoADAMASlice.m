%%%%%%% Fig 10: plot ADAMA (& Litho) Vs slices at a depth
%%%%%%% Fig 11: plot ADAMA velocity slice and velocity gradient at a depth
%%%%%%% Siyu Xue -- May 10. 2023
clear

%% (1) Read in data
load('/Users/sxue3/Documents/BayMap_Figures/Data/Litho1D.mat')  % load Litho 1D models
load('/Users/sxue3/Documents/BayMap_Figures/Data/ADAMAinvert/CADAMA_fixed.mat')  % load ADAMA inversion results
load('/Users/sxue3/Documents/BayMap_Figures/Data/AfrDownSampleDots.mat')  % load data locations
% load the coast
afcoast = load('./Data/AfricaCoast.mat');
aflon = wrapTo180(afcoast.XY(:, 1));
aflat = afcoast.XY(:,2);
coastline = polyshape(aflon, aflat);
% load African country boundaries
S1 = load('./Data/AfrCountry.mat').S1;
madgline = polyshape((S1(448).X).',(S1(448).Y).');  % load Madgascar's Big Island

% load cratons
load('/Users/sxue3/Documents/BayMap_Figures/Data/AfricanCratons.mat');

%% Assemble Litho slice
iD = 51;     % depth = 500*(iD-1) [m]

xgrid = -20:0.1:55;    % lon
ygrid = -35:0.1:40;    % lat
[xq, yq] = meshgrid(xgrid, ygrid);

% Assemble Litho slice
vel2D_litho = zeros(size(dslon));

for ip = 1:size(dslon,1) 
    Lithovs = Litho1D{ip}.pvs;
    vel2D_litho(ip) = Lithovs(iD);
end

F = scatteredInterpolant(dslon,dslat,vel2D_litho, 'linear', 'none');
Vq_litho = F(xq,yq);

for ix = 1:size(Vq_litho, 1)
    for iy = 1:size(Vq_litho, 1)
        if isinterior(coastline, [xgrid(ix), ygrid(iy)]) ||  isinterior(madgline, [xgrid(ix), ygrid(iy)])
            continue
        else
            Vq_litho(iy, ix) = NaN; 
        end
    end
end

% Assemble ADAMA slice

vel2D_ADAMA = zeros(size(dslon));

for ip = 1:size(dslon,1) 
    if ~isempty(ADAMA1D{ip})  % use ADAMA if there is an update
        ADAMAvs = ADAMA1D{ip}.vsv_update;
        fupdate = find(ADAMAvs(:,1) ~= 0);
        vel2D_ADAMA(ip) = ADAMAvs(fupdate(end),iD);
    else
        Lithovs = Litho1D{ip}.pvs;
        vel2D_ADAMA(ip) = Lithovs(iD);
    end
end

F2 = scatteredInterpolant(dslon,dslat,vel2D_ADAMA, 'linear', 'none');
Vq_ADAMA = F2(xq,yq);

for ix = 1:size(Vq_ADAMA, 1)
    for iy = 1:size(Vq_ADAMA, 1)
        if isinterior(coastline, [xgrid(ix), ygrid(iy)]) ||  isinterior(madgline, [xgrid(ix), ygrid(iy)])
            continue
        else
            Vq_ADAMA(iy, ix) = NaN; 
        end
    end
end


%% Plot both maps
t1 = tiledlayout(1,2);
t1.TileSpacing = 'compact';
t1.Padding = 'compact';
clim = [3.5, 4.2];  % colorbar limit

% create a colormap to show the difference
Red = [1,0,0];
White = [1,1,1];
Blue = [0.0745, 0.6235, 1];
RWBcolormap = [[ones(1,128), linspace(1,0.0745,128)]; [linspace(0, 1, 128), linspace(1,0.6235,128)]; [linspace(0, 1, 128), ones(1,128)]].';


ax1 = nexttile;
% Plot the ADAMA slice
pcolor(ax1, xq,yq, Vq_ADAMA/1000)
shading flat
c1 = colorbar;
c1.Title.String = 'Velocity(km/s)';
caxis(clim)
colormap(ax1, RWBcolormap)

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
% title(['ADAMA at depth = ', num2str(500*(iD-1)), ' m'])
set(gca,'fontsize', 14, 'Color', 'w') 

% the second plot
ax2 = nexttile;
% % Plot the Litho slice
% pcolor(xq,yq, Vq_litho/1000)
% shading flat
% c2 = colorbar;
% c2.Title.String = 'Velocity(km/s)';
% caxis(clim)
% colormap(RWBcolormap)
% title(['Litho1.0 at depth = ', num2str(500*(iD-1)), ' m'])

% Plot the ADAMA velocity gradient
load('./Utils/tolu_stdmap2.mat')    % load the colormap
[px, py] = gradient(Vq_ADAMA);
sqGrad = sqrt(px.*px + py.*py);
pcolor(ax2, xq,yq, sqGrad*10/1000)      % [(m/s)/degree]
shading flat
c2 = colorbar('southoutside');
colormap(ax2, ved_cmapa)
caxis([0 .6])
% c2.Title.String = '(km/degree)';
% title('Velocity Gradient')
xlabel('(km/degree)')

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

% set the size of the plot
x0=10;
y0=10;
fwidth=1100;
fheight=700;
set(gcf,'position',[x0,y0,fwidth,fheight], 'Color', 'w');
set(gca,'fontsize', 14) 

% save the figure
% figpath = '/Users/sxue3/Documents/BayMap_Figures/fig/fig11_VelGradient1.pdf';
% set(gcf,'Units','Inches');
% pos = get(gcf,'Position');
% set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
% print(gcf,figpath,'-painters', '-dpdf','-r0');
