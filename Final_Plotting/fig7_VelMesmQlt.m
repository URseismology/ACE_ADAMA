%%%%% Fig 7 (the green-red dots plot): 
%%%%% plot the quality of velocity measurements at all downsampled points
%%%%% Siyu Xue -- Apr 10. 2023

addpath('/Users/sxue3/Documents/BayMap_Figures/Utils')
clear

%% Create sample points in Africa and in cratons
load('/Users/sxue3/Documents/BayMap_Figures/Data/RayleighPData/plotpara.mat');
longrid = longrid(1,:);
latgrid = latgrid(:,1);
% load the coast
afcoast = load('/Users/sxue3/Documents/BayMap_Figures/Data/AfricaCoast.mat');
aflon = wrapTo180(afcoast.XY(:, 1));
aflat = afcoast.XY(:,2);
coastline = polyshape(aflon, aflat);
% load Madgascar's Big Island
S1 = load('/Users/sxue3/Documents/BayMap_Figures/Data/AfrCountry.mat').S1;
madgline = polyshape((S1(448).X).',(S1(448).Y).');  

% load cratons
load('/Users/sxue3/Documents/BayMap_Figures/Data/AfricanCratons.mat');

%%% Old downsample pre-processing code
% % Load all points generated using equal distance
% dsdots = load('/Users/sxue3/Documents/BayMap_Figures/Data/downsample_dots.mat');
% 
% inafr = isinterior(coastline, dsdots.allon, dsdots.allat) + isinterior(madgline, dsdots.allon, dsdots.allat);
% inafr = logical(inafr);
% 
% dslat = dsdots.allat(inafr);
% dslon = dsdots.allon(inafr);
% dslat(1514) = [];  % 1514 is out of the map range!
% dslon(1514) = [];

% % Get the sample points in map
% outdot = [];
% afrmask = zeros(2381, 2);
% for i = 1:length(dslat)
%     xi = find(latgrid <= dslat(i), 1, 'last');
%     yi = find(longrid <= dslon(i), 1, 'last');
% 
%     afrmask(i, :) = [xi, yi];
% end

%%% Now, we can load processed downsample locations:
dsdots = load('/Users/sxue3/Documents/BayMap_Figures/Data/AfrDownSampleDots.mat');
dslat = dsdots.dslat;
dslon = dsdots.dslon;
afrmask = dsdots.afrmask;

%% Read in the velocity reference (PREM)

Rayref = readtable('/Users/sxue3/Documents/BayMap_Figures/Data/SDISPR.ASC', 'FileType','fixedwidth');
Loveref = readtable('/Users/sxue3/Documents/BayMap_Figures/Data/SDISPL.ASC', 'FileType','fixedwidth');
RayGref = readtable('/Users/sxue3/Documents/BayMap_Figures/Data/SREGN.ASC', 'FileType','fixedwidth');

% periods
Prdref = [Rayref{78,3} Rayref{86,3} Rayref{97,3} Rayref{103,3} Rayref{108,3} ...
    Rayref{112,3} Rayref{116,3} Rayref{119,3} Rayref{121,3} Rayref{122,3} Rayref{123,3}];

% Rayleigh phase
Rvelref = [Rayref{78,5} Rayref{86,5} Rayref{97,5} Rayref{103,5} Rayref{108,5} ...
    Rayref{112,5} Rayref{116,5} Rayref{119,5} Rayref{121,5} Rayref{122,5} Rayref{123,5}];

% Love phase
Lvelref = [Loveref{86,5} Loveref{97,5} Loveref{103,5} Loveref{108,5} ...
    Loveref{112,5} Loveref{116,5} Loveref{121,5} Loveref{122,5} Loveref{123,5}];

% Rayleigh group
RGvelref = [RayGref{78,5} RayGref{86,5} RayGref{97,5} RayGref{103,5} RayGref{108,5} ...
    RayGref{112,5} RayGref{116,5} RayGref{119,5} RayGref{121,5} RayGref{122,5} RayGref{123,5}];

% Love group
LGvelref = [3.68, 3.91, 4.07, 4.18]; % for now, we only have data at 25s, 30s, 35s, and 40s

%% Get the downsampled velocity
type = 'G'; % 'P' or 'G'
% Periods = {'R5', 'R6', 'R8', 'R10', 'R12', 'R15', 'R20', 'R25', 'R30', 'R35', 'R40'};
% Periods = {'L6', 'L8', 'L10', 'L12', 'L15', 'L20', 'L30', 'L35', 'L40'};
Periods = {'L25', 'L30', 'L35', 'L40'};
% velref = Rvelref;
% velref = RGvelref;
% velref = Lvelref;
velref = LGvelref;

% set the standards
SDcut = 3;      % more than this is high SD
BIAScut = 2;    % more than this is biased

for ip = 1:length(Periods)
    period = Periods{ip};
    
%     pavgpath = ['/Users/sxue3/Documents/BayMap_Figures/Data/RayAvgMap/', period, '_', type, '_maps.mat'];
    pavgpath = ['/Users/sxue3/Documents/BayMap_Figures/Data/LoveAvgMap/', period, '_', type, '_maps.mat'];

    pavg = load(pavgpath).avgmap;
    psd = load(pavgpath).sdevmap;

    ipavg = zeros(length(dslat), 1);
    ipsd = zeros(length(dslat), 1);
    for id = 1:length(dslat)
        ipavg(id) = pavg(afrmask(id, 1), afrmask(id, 2));
        ipsd(id) = psd(afrmask(id, 1), afrmask(id, 2));
    end 
    pavgs(:,ip) = ipavg;
    psds(:,ip) = ipsd;
end

%% Get the locations of good points and bad points

hold on
plot(coastline, 'EdgeColor', 'k', 'LineWidth', 2, 'FaceAlpha', 0);  
% Plot Madgascar's Big Island
plot(madgline.Vertices(1:50:end, 1), madgline.Vertices(1:50:end, 2),'k', 'LineWidth', 2);

%%%% Loop through all down sampled points and pick the bad curves at 6s
bvhsd = [];
bvlsd = [];
gvhsd = [];
gvlsd = [];
totpts = length(dslat);

% psd_threshold = mean(psds,1);
psd_threshold = ones(length(Periods), 1)*0.4;
dp_record = zeros(length(dslat),1);

for id = 1:length(dslat)

    if sum(psds(id, :) >= psd_threshold) > SDcut
        highSD = 1;
    else
        highSD = 0;
    end

    velout = sum(pavgs(id,:) < velref*0.6) + sum(pavgs(id,:) > velref*1.4);
    if velout > BIAScut  % then the result is biased
        if ~highSD  % if low SD  (red)
            p1 = plot(dslon(id), dslat(id), '.', 'MarkerSize',15, 'Color',[1 0 0]); 
            bvlsd = [bvlsd id];
            dp_record(id) = 3;
        else  % if high SD  (dark red)
            p2 = plot(dslon(id), dslat(id), '.', 'MarkerSize',15, 'Color',[0.6353 0.0784 0.1843]); 
            bvhsd = [bvhsd id];
            dp_record(id) = 4;
        end
    else
        if ~highSD  % if low SD  (green)
            p3 = plot(dslon(id), dslat(id), '.', 'MarkerSize',15, 'Color',[0.1059 0.9804 0.1059]);
            gvlsd = [gvlsd id];
            dp_record(id) = 1;
        else  % if high SD  (blue)
            p4 = plot(dslon(id), dslat(id), '.', 'MarkerSize',15, 'Color',[0.3020 0.7451 0.9333]);
            gvhsd = [gvhsd  id];
            dp_record(id) = 2;
        end
    end
end

% get the % of the green dots in the map
% length(gvlsd)*100 / (length(bvhsd) + length(bvlsd) + length(gvlsd) + length(gvhsd))

% Plot African cratons
plot(AC.Congo(:,2), AC.Congo(:,1), 'color', 'k', 'lineWidth', 2);
plot(AC.Kala(:,2), AC.Kala(:,1), 'color', 'k', 'lineWidth', 2);
plot(AC.Sahara(:,2), AC.Sahara(:,1), 'color', 'k', 'lineWidth', 2);
plot(AC.Tanza(:,2), AC.Tanza(:,1), 'color', 'k', 'lineWidth', 2);
plot(AC.West(:,2), AC.West(:,1), 'color', 'k', 'lineWidth', 2);
plot(AC.Zim(:,2), AC.Zim(:,1), 'color', 'k', 'lineWidth', 2);

% % Only plot in Rayleigh Phase map
% % add the transec in Fig. 9
% xs = 0:0.1:50;
% ys1 = 24/25 * xs - 21.5;
% plot(xs, ys1, '--', 'Color', [0.502 0.502 0.502], 'LineWidth',3)
% % add the examples in Fig. 8
% plot(dslon(7), dslat(7),'s','MarkerSize',10,'MarkerEdgeColor',[1 0.7725 0.2392],'lineWidth', 2)
% plot(dslon(6), dslat(6),'s','MarkerSize',10,'MarkerEdgeColor',[1 0.7725 0.2392],'lineWidth', 2)
% plot(dslon(28), dslat(28),'s','MarkerSize',10,'MarkerEdgeColor',[1 0.7725 0.2392],'lineWidth', 2)
% plot(dslon(29), dslat(29),'s','MarkerSize',10,'MarkerEdgeColor',[1 0.7725 0.2392],'lineWidth', 2)

% plot settings
minlon = min(longrid);
maxlon = max(longrid);
minlat = min(latgrid);
maxlat = max(latgrid);
xlim([minlon maxlon]);
ylim([minlat maxlat]);
axesm('mercator','maplatlimit',[minlat maxlat],'maplonlimit',[minlon maxlon]);
xticks(10*ceil(minlon/10):10:10*floor(maxlon/10));
yticks(10*ceil(minlat/10):10:10*floor(maxlat/10));

legend([p3 p4 p1 p2],{'High precision (unbiased)', 'Low precision (unbiased)'...
    'Biased (high precision)', 'Biased (low precision)'}, 'Location','bestoutside')
% title('Rayleigh - Phase')

% set the size of the plot
x0=10;
y0=10;
fwidth=800;
fheight=500;
set(gcf,'position',[x0,y0,fwidth,fheight]);

% % save the figure
% figpath = '/Users/sxue3/Documents/BayMap_Figures/fig/fig7_RayP_mapqlt.pdf';
% set(gcf,'Units','Inches');
% pos = get(gcf,'Position');
% set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
% print(gcf,figpath,'-painters', '-dpdf','-r0');


