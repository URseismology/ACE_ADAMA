%%%%% Fig Suppliment
%%%%% Show the location where Litho model is updated

clear

%% Load the data
load('/Users/sxue3/Documents/BayMap_Figures/Data/ADAMAinvert/ADAMA_1D_fixed1p7.mat')  % load ADAMA inversion results
load('/Users/sxue3/Documents/BayMap_Figures/Data/AfrDownSampleDots.mat')  % load data locations
% load the coast
afcoast = load('./Data/GeoData/AfricaCoast.mat');
aflon = wrapTo180(afcoast.XY(:, 1));
aflat = afcoast.XY(:,2);
coastline = polyshape(aflon, aflat);
% load African country boundaries
S1 = load('./Data/GeoData/AfrCountry.mat').S1;
madgline = polyshape((S1(448).X).',(S1(448).Y).');  % load Madgascar's Big Island

% load cratons
load('/Users/sxue3/Documents/BayMap_Figures/Data/GeoData/AfricanCratons.mat');


%% Plot

hold on
totalUpdate = 0;
for ip = 1:size(dslon,1) 
    if ~isempty(ADAMA1D{ip})  % use ADAMA if there is an update
        plot(dslon(ip), dslat(ip), 'r.', 'MarkerSize', 10)
        totalUpdate = totalUpdate + 1;
    end
end


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

% % save the figure
% figpath = '/Users/sxue3/Documents/BayMap_Figures/fig/figSX_LithoUpdateMap.pdf';
% set(gcf,'Units','Inches');
% pos = get(gcf,'Position');
% set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
% print(gcf,figpath,'-painters', '-dpdf','-r0');


