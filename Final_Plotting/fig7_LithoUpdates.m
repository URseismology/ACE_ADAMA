%%%% Plot the Litho updates at sampled locations

%% Read in Data

% read Litho model
load('/Users/sxue3/Documents/BayMap_Figures/Data/Litho1D.mat');
load('/Users/sxue3/Documents/BayMap_Figures/Data/Litho_dwsample.mat');
% read ADAMA Litho updates
load('/Users/sxue3/Documents/BayMap_Figures/Data/ADAMAinvert/ADAMA_1D_fixed1p7.mat');

% read sample coordinates
dsdots = load('/Users/sxue3/Documents/BayMap_Figures/Data/AfrDownSampleDots.mat');
dslat = dsdots.dslat;
dslon = dsdots.dslon;
afrmask = dsdots.afrmask;

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
lonZim = AC.Zim(:,2);% load the coast
afcoast = load('/Users/sxue3/Documents/BayMap_Figures/Data/AfricaCoast.mat');
aflon = wrapTo180(afcoast.XY(:, 1));
aflat = afcoast.XY(:,2);
coastline = polyshape(aflon, aflat);
% load Madgascar's Big Island
S1 = load('/Users/sxue3/Documents/BayMap_Figures/Data/AfrCountry.mat').S1;
madgline = polyshape((S1(448).X).',(S1(448).Y).');  

% model we are interesed in
im1 = 7;
im2 = 8;
im3 = 19;
im4 = 20;

% Load ADAMA velocity measurement
type = 'P'; % 'P' or 'G'
Periods = {'R40', 'R35', 'R30', 'R25', 'R20', 'R15', 'R12', 'R10', 'R8', 'R6', 'R5'};
% Periods = {'L5', 'L6', 'L8', 'L10', 'L12', 'L15', 'L20', 'L25', 'L30', 'L35', 'L40'};

for ip = 1:length(Periods)
    period = Periods{ip};
    
    pavgpath = ['/Users/sxue3/Documents/BayMap_Figures/Data/RayAvgMap/', period, '_', type, '_maps.mat'];
%     pavgpath = ['/Users/sxue3/Documents/BayMap_Figures/Data/LoveAvgMap/', period, '_', type, '_maps.mat'];
    pavg = load(pavgpath).avgmap;
    psd = load(pavgpath).sdevmap;

    % get the vel. and sd. at a point
    pavg1(ip) = pavg(afrmask(im1, 1), afrmask(im1, 2));
    psd1(ip) = psd(afrmask(im1, 1), afrmask(im1, 2));
    pavg2(ip) = pavg(afrmask(im2, 1), afrmask(im2, 2));
    psd2(ip) = psd(afrmask(im2, 1), afrmask(im2, 2));
    pavg3(ip) = pavg(afrmask(im3, 1), afrmask(im3, 2));
    psd3(ip) = psd(afrmask(im3, 1), afrmask(im3, 2));
    pavg4(ip) = pavg(afrmask(im4, 1), afrmask(im4, 2));
    psd4(ip) = psd(afrmask(im4, 1), afrmask(im4, 2));
end

%% Plot everything
set(0,'defaultfigurecolor',[1 1 1])

%%%%%%%%%%%%%%%%%% PLOT dispersion curves %%%%%%%%%%%%%%%%%%
ax11 = axes('position',[0.1, 0.79, 0.2, 0.18]);
ax12 = axes('position',[0.1, 0.56, 0.2, 0.18]);
ax13 = axes('position',[0.1, 0.33, 0.2, 0.18]);
ax14 = axes('position',[0.1, 0.1, 0.2, 0.18]);

hold(ax11,'on')
p1 = plot(ax11,LithoRay(im1, :), 'k');   % Litho rayleigh velocity
p2 = plot(ax11,flip(ADAMA1D{im1}.U), 'r');     % the final ADAMA update
p3 = errorbar(ax11,pavg1.*1000, psd1.*1000, 'b');    % ADAMA measurements w. uncert.
xticks(ax11,[1 2 3 4 5 6 7 8 9 10 11])
xticklabels(ax11,[40 35 30 25 20 15 12 10 8 6 5]);
% xlabel('Period (s)')
ylabel(ax11,'Velocity (m/s)')
legend(ax11,'Litho', 'final', 'ADAMA', 'Location','southwest')
hold(ax11,'off')

hold(ax12,'on')
p1 = plot(ax12,LithoRay(im2, :), 'k');   % Litho rayleigh velocity
p2 = plot(ax12,flip(ADAMA1D{im2}.U), 'r');     % the final ADAMA update
p3 = errorbar(ax12,pavg2.*1000, psd2.*1000, 'b');    % ADAMA measurements w. uncert.
xticks(ax12,[1 2 3 4 5 6 7 8 9 10 11])
xticklabels(ax12,[40 35 30 25 20 15 12 10 8 6 5]);
% xlabel('Period (s)')
ylabel(ax12,'Velocity (m/s)')
hold(ax12,'off')

hold(ax13,'on')
p1 = plot(ax13,LithoRay(im3, :), 'k');   % Litho rayleigh velocity
p2 = plot(ax13,flip(ADAMA1D{im3}.U), 'r');     % the final ADAMA update
p3 = errorbar(ax13,pavg3.*1000, psd3.*1000, 'b');    % ADAMA measurements w. uncert.
xticks(ax13,[1 2 3 4 5 6 7 8 9 10 11])
xticklabels(ax13,[40 35 30 25 20 15 12 10 8 6 5]);
% xlabel('Period (s)')
ylabel(ax13,'Velocity (m/s)')
hold(ax13,'off')

hold(ax14,'on')
p1 = plot(ax14,LithoRay(im4, :), 'k');   % Litho rayleigh velocity
p2 = plot(ax14,flip(ADAMA1D{im4}.U), 'r');     % the final ADAMA update
p3 = errorbar(ax14,pavg4.*1000, psd4.*1000, 'b');    % ADAMA measurements w. uncert.
xticks(ax14,[1 2 3 4 5 6 7 8 9 10 11])
xticklabels(ax14,[40 35 30 25 20 15 12 10 8 6 5]);
xlabel(ax14,'Period (s)')
ylabel(ax14,'Velocity (m/s)')
hold(ax14,'off')

%%%%%%%%%%%%%%%%%% PLOT models and model updates %%%%%%%%%%%%%%%%%%
ax21 = axes('position',[0.35, 0.55, 0.1, 0.4]);
ax22 = axes('position',[0.5, 0.55, 0.1, 0.4]);
ax23 = axes('position',[0.35, 0.1, 0.1, 0.4]);
ax24 = axes('position',[0.5, 0.1, 0.1, 0.4]);

% convert Litho to the same format
vSGrid = Litho1D{im1}.vs;
thckGrid = Litho1D{im1}.thick;
% in case we get fluid layer, delete it
if vSGrid(1) == 0
    thckGrid(2) = thckGrid(2) + thckGrid(1);
    thckGrid(1) = [];
end
maxmodel = 600;
thckGrid = thckGrid(thckGrid ~= 0);
vSGrid = vSGrid(vSGrid ~= 0);
% compute the depth of each layer
depthGrid = zeros(length(thckGrid), 1);
for id = flip(1:length(thckGrid))
    depthGrid(id) = sum(thckGrid(1:id));
end
depthGrid = round(depthGrid./500);
if sum(depthGrid) > 600
    % ignore the model deeper than 300m
    depthGrid(end) = 600 - sum(depthGrid(1:end-1));
end
% create the input for Raylee code
pvs = ones(maxmodel, 1);
layerst = 1;
for i = 1:length(depthGrid)
    pvs(layerst: depthGrid(i)) = pvs(layerst: depthGrid(i)) .* vSGrid(i);
    layerst = 1 + depthGrid(i);
end
% fill to the max. depth
pvs(layerst: end) = pvs(layerst: end) .* vSGrid(end);
% get the updated Litho model
ADAMAvs = ADAMA1D{im1}.vsv_update;
fupdate = find(ADAMAvs(:,1) ~= 0);
% plot two models
hold(ax21,'on')
plot(ax21, ADAMAvs(fupdate(end),:), 'r')
plot(ax21, pvs, 'k')
xlim(ax21, [1 120])
xticks(ax21, [0 20 40 60 80 100 120])
xticklabels(ax21, {'0','10','20','30','40','50','60'})  % grid size is 500m
camroll(ax21, -90)
xlabel(ax21, 'Depth (km)')
ylabel(ax21, 'Velocity (m/s)')
ylim(ax21, [0 5000])
hold(ax21,'off')

% convert Litho to the same format
vSGrid = Litho1D{im2}.vs;
thckGrid = Litho1D{im2}.thick;
% in case we get fluid layer, delete it
if vSGrid(1) == 0
    thckGrid(2) = thckGrid(2) + thckGrid(1);
    thckGrid(1) = [];
end
maxmodel = 600;
thckGrid = thckGrid(thckGrid ~= 0);
vSGrid = vSGrid(vSGrid ~= 0);
% compute the depth of each layer
depthGrid = zeros(length(thckGrid), 1);
for id = flip(1:length(thckGrid))
    depthGrid(id) = sum(thckGrid(1:id));
end
depthGrid = round(depthGrid./500);
if sum(depthGrid) > 600
    % ignore the model deeper than 300m
    depthGrid(end) = 600 - sum(depthGrid(1:end-1));
end
% create the input for Raylee code
pvs = ones(maxmodel, 1);
layerst = 1;
for i = 1:length(depthGrid)
    pvs(layerst: depthGrid(i)) = pvs(layerst: depthGrid(i)) .* vSGrid(i);
    layerst = 1 + depthGrid(i);
end
% fill to the max. depth
pvs(layerst: end) = pvs(layerst: end) .* vSGrid(end);
% get the updated Litho model
ADAMAvs = ADAMA1D{im2}.vsv_update;
fupdate = find(ADAMAvs(:,1) ~= 0);
% plot two models
hold(ax22,'on')
plot(ax22, pvs, 'k')
plot(ax22, ADAMAvs(fupdate(end),:), 'r')
xlim(ax22, [1 120])
xticks(ax22, [0 20 40 60 80 100 120])
xticklabels(ax22, {'0','10','20','30','40','50','60'})  % grid size is 500m
camroll(ax22, -90)
% xlabel(ax22, 'Depth (km)')
ylabel(ax22, 'Velocity (m/s)')
ylim(ax22, [0 5000])
hold(ax22,'off')

% convert Litho to the same format
vSGrid = Litho1D{im3}.vs;
thckGrid = Litho1D{im3}.thick;
% in case we get fluid layer, delete it
if vSGrid(1) == 0
    thckGrid(2) = thckGrid(2) + thckGrid(1);
    thckGrid(1) = [];
end
maxmodel = 600;
thckGrid = thckGrid(thckGrid ~= 0);
vSGrid = vSGrid(vSGrid ~= 0);
% compute the depth of each layer
depthGrid = zeros(length(thckGrid), 1);
for id = flip(1:length(thckGrid))
    depthGrid(id) = sum(thckGrid(1:id));
end
depthGrid = round(depthGrid./500);
if sum(depthGrid) > 600
    % ignore the model deeper than 300m
    depthGrid(end) = 600 - sum(depthGrid(1:end-1));
end
% create the input for Raylee code
pvs = ones(maxmodel, 1);
layerst = 1;
for i = 1:length(depthGrid)
    pvs(layerst: depthGrid(i)) = pvs(layerst: depthGrid(i)) .* vSGrid(i);
    layerst = 1 + depthGrid(i);
end
% fill to the max. depth
pvs(layerst: end) = pvs(layerst: end) .* vSGrid(end);
% get the updated Litho model
ADAMAvs = ADAMA1D{im3}.vsv_update;
fupdate = find(ADAMAvs(:,1) ~= 0);
% plot two models
hold(ax23,'on')
plot(ax23, pvs, 'k')
plot(ax23, ADAMAvs(fupdate(end),:), 'r')
xlim(ax23, [1 120])
xticks(ax23, [0 20 40 60 80 100 120])
xticklabels(ax23, {'0','10','20','30','40','50','60'})  % grid size is 500m
camroll(ax23, -90)
xlabel(ax23, 'Depth (km)')
% ylabel(ax23, 'Velocity (m/s)')
ylim(ax23, [0 5000])
hold(ax23,'off')

% convert Litho to the same format
vSGrid = Litho1D{im4}.vs;
thckGrid = Litho1D{im4}.thick;
% in case we get fluid layer, delete it
if vSGrid(1) == 0
    thckGrid(2) = thckGrid(2) + thckGrid(1);
    thckGrid(1) = [];
end
maxmodel = 600;
thckGrid = thckGrid(thckGrid ~= 0);
vSGrid = vSGrid(vSGrid ~= 0);
% compute the depth of each layer
depthGrid = zeros(length(thckGrid), 1);
for id = flip(1:length(thckGrid))
    depthGrid(id) = sum(thckGrid(1:id));
end
depthGrid = round(depthGrid./500);
if sum(depthGrid) > 600
    % ignore the model deeper than 300 km
    depthGrid(end) = 600 - sum(depthGrid(1:end-1));
end
% create the input for Raylee code
pvs = ones(maxmodel, 1);
layerst = 1;
for i = 1:length(depthGrid)
    pvs(layerst: depthGrid(i)) = pvs(layerst: depthGrid(i)) .* vSGrid(i);
    layerst = 1 + depthGrid(i);
end
% fill to the max. depth
pvs(layerst: end) = pvs(layerst: end) .* vSGrid(end);
% get the updated Litho model
ADAMAvs = ADAMA1D{im4}.vsv_update;
fupdate = find(ADAMAvs(:,1) ~= 0);
% plot two models
hold(ax24,'on')
plot(ax24, pvs, 'k')
plot(ax24, ADAMAvs(fupdate(end),:), 'r')
xlim(ax24, [1 120])
xticks(ax24, [0 20 40 60 80 100 120])
xticklabels(ax24, {'0','10','20','30','40','50','60'})  % grid size is 500m
camroll(ax24, -90)
% xlabel(ax23, 'Depth (km)')
% ylabel(ax23, 'Velocity (m/s)')
ylim(ax24, [0 5000])
hold(ax24,'off')


%%%%%%%%%%%%%%%%%% PLOT African map with locations %%%%%%%%%%%%%%%%%%
ax3 = axes('position',[0.64, 0.2, 0.35, 0.7]);

hold on
plot(ax3, coastline, 'EdgeColor', 'k', 'LineWidth', 2, 'FaceAlpha', 0);  
% Plot Madgascar's Big Island
plot(ax3, madgline.Vertices(1:50:end, 1), madgline.Vertices(1:50:end, 2),'k', 'LineWidth', 2);
% Plot African cratons
plot(ax3, lonCongo, latCongo, 'color', 'k', 'lineWidth', 2);
plot(ax3, lonKala, latKala, 'color', 'k', 'lineWidth', 2);
plot(ax3, lonSahara, latSahara, 'color', 'k', 'lineWidth', 2);
plot(ax3, lonTanza, latTanza, 'color', 'k', 'lineWidth', 2);
plot(ax3, lonWest, latWest, 'color', 'k', 'lineWidth', 2);
plot(ax3, lonZim, latZim, 'color', 'k', 'lineWidth', 2);

% plot the dots
d1 = plot(ax3, dslon(im1), dslat(im1), '.', 'MarkerSize',20, 'DisplayName','sample 1');
d2 = plot(ax3, dslon(im2), dslat(im2), '.', 'MarkerSize',20, 'DisplayName','sample 2');
d3 = plot(ax3, dslon(im3), dslat(im3), '.', 'MarkerSize',20, 'DisplayName','sample 3');
d4 = plot(ax3, dslon(im4), dslat(im4), '.', 'MarkerSize',20, 'DisplayName','sample 4');

legend([d1 d2 d3 d4])

% plot settings
minlon = -20;
maxlon = 55;
minlat = -35;
maxlat = 40;
xlim([minlon maxlon]);
ylim([minlat maxlat]);
axesm('mercator','maplatlimit',[minlat maxlat],'maplonlimit',[minlon maxlon]);
xticks(10*ceil(minlon/10):10:10*floor(maxlon/10));
yticks(10*ceil(minlat/10):10:10*floor(maxlat/10));
hold off


% set the size of the plot
x0=10;
y0=10;
fwidth=1200;
fheight=600;
set(gcf,'position',[x0,y0,fwidth,fheight]);