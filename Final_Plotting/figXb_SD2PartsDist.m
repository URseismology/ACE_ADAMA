%%%% Plot the boxplot to show the distribution of velocity SD and percentage error in Africa
addpath('/Users/sxue3/Documents/BayMap_Figures/Utils')
%% Read in the data
% read in the SD in Africa
% load the coast
load('/Users/sxue3/Documents/BayMap_Figures/Data/fig2P/plotpara.mat');
afcoast = load('/Users/sxue3/Documents/BayMap_Figures/Data/AfricaCoast.mat');
aflon = wrapTo180(afcoast.XY(:, 1));
aflat = afcoast.XY(:,2);
coastline = polyshape(aflon, aflat);
% load Madgascar's Big Island
S1 = load('/Users/sxue3/Documents/BayMap_Figures/Data/AfrCountry.mat').S1;
madgline = polyshape((S1(448).X).',(S1(448).Y).');  

afrmask = zeros(297, 324);
mapsize = size(afrmask);
for xi = 1:mapsize(2)
    loni = longrid(1,xi);
    for yi = 1:mapsize(1)
        lati = latgrid(yi,1);

        if isinterior(coastline, loni, lati) == 1 || isinterior(madgline, loni, lati) == 1
            afrmask(yi, xi) = 1;
        end
    end
end

afrmask = logical(afrmask);

PSD40 = load('/Users/sxue3/Documents/BayMap_Figures/Data/fig4maps/R40_P.mat').sdevmap;
PSD40 = PSD40(afrmask);

GSD40 = load('/Users/sxue3/Documents/BayMap_Figures/Data/fig4maps/R40_G.mat').sdevmap;
GSD40 = GSD40(afrmask);

PSD20 = load('/Users/sxue3/Documents/BayMap_Figures/Data/fig4maps/R20_P.mat').sdevmap;
PSD20 = PSD20(afrmask);

GSD20 = load('/Users/sxue3/Documents/BayMap_Figures/Data/fig4maps/R20_G.mat').sdevmap;
GSD20 = GSD20(afrmask);

PSD15 = load('/Users/sxue3/Documents/BayMap_Figures/Data/fig4maps/R15_P.mat').sdevmap;
PSD15 = PSD15(afrmask);

GSD15 = load('/Users/sxue3/Documents/BayMap_Figures/Data/fig4maps/R15_G.mat').sdevmap;
GSD15 = GSD15(afrmask);

PSD10 = load('/Users/sxue3/Documents/BayMap_Figures/Data/fig4maps/R10_P.mat').sdevmap;
PSD10 = PSD10(afrmask);

GSD10 = load('/Users/sxue3/Documents/BayMap_Figures/Data/fig4maps/R10_G.mat').sdevmap;
GSD10 = GSD10(afrmask);

PSD5 = load('/Users/sxue3/Documents/BayMap_Figures/Data/fig4maps/R5_P.mat').sdevmap;
PSD5 = PSD5(afrmask);

GSD5 = load('/Users/sxue3/Documents/BayMap_Figures/Data/fig4maps/R5_G.mat').sdevmap;
GSD5 = GSD5(afrmask);

%% get the top and bottom half of the SD

PSD40 = sort(PSD40);
bPSD40 = PSD40(1:floor(length(PSD40)/2));
tPSD40 = PSD40(ceil(length(PSD40)/2):end);

PSD20 = sort(PSD20);
bPSD20 = PSD20(1:floor(length(PSD20)/2));
tPSD20 = PSD20(ceil(length(PSD20)/2):end);

PSD15 = sort(PSD15);
bPSD15 = PSD15(1:floor(length(PSD15)/2));
tPSD15 = PSD15(ceil(length(PSD15)/2):end);

PSD10 = sort(PSD10);
bPSD10 = PSD10(1:floor(length(PSD10)/2));
tPSD10 = PSD10(ceil(length(PSD10)/2):end);

PSD5 = sort(PSD5);
bPSD5 = PSD5(1:floor(length(PSD5)/2));
tPSD5 = PSD5(ceil(length(PSD5)/2):end);

GSD40 = sort(GSD40);
bGSD40 = GSD40(1:floor(length(GSD40)/2));
tGSD40 = GSD40(ceil(length(GSD40)/2):end);

GSD20 = sort(GSD20);
bGSD20 = GSD20(1:floor(length(GSD20)/2));
tGSD20 = GSD20(ceil(length(GSD20)/2):end);

GSD15 = sort(GSD15);
bGSD15 = GSD15(1:floor(length(GSD15)/2));
tGSD15 = GSD15(ceil(length(GSD15)/2):end);

GSD10 = sort(GSD10);
bGSD10 = GSD10(1:floor(length(GSD10)/2));
tGSD10 = GSD10(ceil(length(GSD10)/2):end);

GSD5 = sort(GSD5);
bGSD5 = GSD5(1:floor(length(GSD5)/2));
tGSD5 = GSD5(ceil(length(GSD5)/2):end);
%% initialize some parameters
numgroups = 5;
DeviceColors = {[0.59 0.89 1],[0 0.45 0.74],[1 0.82 0.46], [1 0.41 0.01]};
legendEntries = {'P-low','G-low', 'P-high','G-high'};
% Periods = {'40' '20' '15' '10' '05'};
% Freqs = {'25' '50' '67' '100' '200'};
Periods = {'40' '20' '15' '10' '05'};
Freqs = {'25' '50' '67' '100' '200'};
N = 4;  % #of boxplot in a group
delta = linspace(-.35,.35,5); %// define offsets to distinguish plots
width = eps; %// small width to avoid overlap
legWidth = 2.3; %// make room for legend

grp1data = [bPSD40.' bPSD20.' bPSD15.' bPSD10.' bPSD5.'];
grp1idx = [(zeros(length(bPSD40),1)+1).' (zeros(length(bPSD20),1)+2).' ...
    (zeros(length(bPSD15),1)+3).' (zeros(length(bPSD10),1)+4).' (zeros(length(bPSD5),1)+5).'];

grp2data = [bGSD40.' bGSD20.' bGSD15.' bGSD10.' bGSD5.'];
grp2idx = [(zeros(length(bGSD40),1)+1).' (zeros(length(bGSD20),1)+2).' ...
    (zeros(length(bGSD15),1)+3).' (zeros(length(bGSD10),1)+4).' (zeros(length(bGSD5),1)+5).'];

grp3data = [tPSD40.' tPSD20.' tPSD15.' tPSD10.' tPSD5.'];
grp3idx = [(zeros(length(tPSD40),1)+1).' (zeros(length(tPSD20),1)+2).' ...
    (zeros(length(tPSD15),1)+3).' (zeros(length(tPSD10),1)+4).' (zeros(length(tPSD5),1)+5).'];

grp4data = [tGSD40.' tGSD20.' tGSD15.' tGSD10.' tGSD5.'];
grp4idx = [(zeros(length(tGSD40),1)+1).' (zeros(length(tGSD20),1)+2).' ...
    (zeros(length(tGSD15),1)+3).' (zeros(length(tGSD10),1)+4).' (zeros(length(tGSD5),1)+5).'];



% plot all boxplots
for ii=1:N 
    labels = Periods; 
    hold on

    if ii == 1  %% plotting Phase Cell count
        yyaxis left
        boxplot(grp1data, grp1idx,'Color', DeviceColors{ii}, 'boxstyle','filled', ...
        'position',(1:numgroups)+delta(ii), 'widths',width, 'labels',labels, 'symbol', 'o', ...
        'outliersize', 1.0, 'whisker', 0,  'notch', 'off');
        medians = [nanmedian(bPSD40) nanmedian(bPSD20) nanmedian(bPSD15) nanmedian(bPSD10) nanmedian(bPSD5)];
     
        pp = plot((1:numgroups)+delta(ii), medians, 'ro');
        pp.MarkerFaceColor = DeviceColors{ii};
        pp.MarkerEdgeColor = 'k';
        pp.MarkerSize = 10;

    elseif ii == 2  %% plotting Group Cell count
        yyaxis left
        boxplot(grp2data, grp2idx,'Color', DeviceColors{ii}, 'boxstyle','filled', ...
        'position',(1:numgroups)+delta(ii), 'widths',width, 'labels',labels, 'symbol', 'o', ...
        'outliersize', 1.0, 'whisker', 0,  'notch', 'off');
        medians = [nanmedian(bGSD40) nanmedian(bGSD20) nanmedian(bGSD15) nanmedian(bGSD10) nanmedian(bGSD5)];
    
        pp = plot((1:numgroups)+delta(ii), medians, 'ro');
        pp.MarkerFaceColor = DeviceColors{ii};
        pp.MarkerEdgeColor = 'k';
        pp.MarkerSize = 10;

    elseif ii == 3  %% plotting Phase Velocity
        yyaxis right
        boxplot(grp3data, grp3idx,'Color', DeviceColors{ii}, 'boxstyle','filled', ...
        'position',(1:numgroups)+delta(ii), 'widths',width, 'labels',labels, 'symbol', 'o', ...
        'outliersize', 1.0, 'whisker', 0,  'notch', 'off');
        medians = [nanmedian(tPSD40) nanmedian(tPSD20) nanmedian(tPSD15) nanmedian(tPSD10) nanmedian(tPSD5)];
        
        pp_no = plot((1:numgroups)+delta(ii), medians, 'ro');
        pp_no.MarkerFaceColor = DeviceColors{ii};
        pp_no.MarkerEdgeColor = 'k';
        pp_no.MarkerSize = 10;

    elseif ii == 4  %% plotting Group Velocity
        yyaxis right
        boxplot(grp4data, grp4idx,'Color', DeviceColors{ii}, 'boxstyle','filled', ...
        'position',(1:numgroups)+delta(ii), 'widths',width, 'labels',labels, 'symbol', 'o', ...
        'outliersize', 1.0, 'whisker', 0,  'notch', 'off');
        medians = [nanmedian(tGSD40) nanmedian(tGSD20) nanmedian(tGSD15) nanmedian(tGSD10) nanmedian(tGSD5)];
    
        pp_no = plot((1:numgroups)+delta(ii), medians, 'ro');
        pp_no.MarkerFaceColor = DeviceColors{ii};
        pp_no.MarkerEdgeColor = 'k';
        pp_no.MarkerSize = 10;
    end  
        
end

grid on;
legend(legendEntries, 'Location', 'northwest','fontsize',10);

yyaxis left
ylabel('Standard Deviation (bottom half)','fontsize',15);
ylim([0 2]);
% text(0.4, 1.2,'(d)','fontsize',25);
 
yyaxis right
ylabel('Standard Deviation (top half)','fontsize',15);
ylim([0 2]);

xlabel('Period (s)','fontsize',15);
%set(gca,'XAxisLocation','top','LineWidth',1);  % make the x-axis on top 

% xlabel('Frequency (mHz)','fontsize',15);
set(gca,'LineWidth',1);  

% set the size of the plot
x0=10;
y0=10;
fwidth=600;
fheight=400;
set(gcf,'position',[x0,y0,fwidth,fheight]);

% save the figure
% figpath = '/Users/sxue3/Documents/BayMap_Figures/fig/fig5b_SD2PartsDist.pdf';
% set(gcf,'Units','Inches');
% pos = get(gcf,'Position');
% set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
% print(gcf,figpath,'-painters', '-dpdf','-r0');