%%%% Plot the boxplot to show the distribution of velocity SD and percentage error in Africa
addpath('/Users/sxue3/Documents/BayMap_Figures/Utils')
%% Read in the data

% Read in the percentage error 
err_path = '/Users/sxue3/Documents/BayMap_Figures/Data/map_error/observation_errors_R40P.txt';
A = readtable(err_path,'FileType','text');
PE40= A{:,'Var8'};
err_path = '/Users/sxue3/Documents/BayMap_Figures/Data/map_error/observation_errors_R20P.txt';
A = readtable(err_path,'FileType','text');
PE20= A{:,'Var8'};
err_path = '/Users/sxue3/Documents/BayMap_Figures/Data/map_error/observation_errors_R15P.txt';
A = readtable(err_path,'FileType','text');
PE15= A{:,'Var8'};
err_path = '/Users/sxue3/Documents/BayMap_Figures/Data/map_error/observation_errors_R10P.txt';
A = readtable(err_path,'FileType','text');
PE10= A{:,'Var8'};
err_path = '/Users/sxue3/Documents/BayMap_Figures/Data/map_error/observation_errors_R5P.txt';
A = readtable(err_path,'FileType','text');
PE5= A{:,'Var8'};

err_path = '/Users/sxue3/Documents/BayMap_Figures/Data/map_error/observation_errors_R40G.txt';
A = readtable(err_path,'FileType','text');
GE40= A{:,'Var8'};
err_path = '/Users/sxue3/Documents/BayMap_Figures/Data/map_error/observation_errors_R20G.txt';
A = readtable(err_path,'FileType','text');
GE20= A{:,'Var8'};
err_path = '/Users/sxue3/Documents/BayMap_Figures/Data/map_error/observation_errors_R15G.txt';
A = readtable(err_path,'FileType','text');
GE15= A{:,'Var8'};
err_path = '/Users/sxue3/Documents/BayMap_Figures/Data/map_error/observation_errors_R10G.txt';
A = readtable(err_path,'FileType','text');
GE10= A{:,'Var8'};
err_path = '/Users/sxue3/Documents/BayMap_Figures/Data/map_error/observation_errors_R5G.txt';
A = readtable(err_path,'FileType','text');
GE5= A{:,'Var8'};
clear A

% read in the SD in Africa
%% load the coast
load('/Users/sxue3/Documents/BayMap_Figures/Data/RayleighPData/plotpara.mat');
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


%% initialize some parameters
numgroups = 5;
DeviceColors = {[0.59 0.89 1],[0 0.45 0.74],[1 0.82 0.46], [1 0.41 0.01]};
legendEntries = {'P-SD','G-SD', 'P-Err','G-Err'};
% Periods = {'40' '20' '15' '10' '05'};
% Freqs = {'25' '50' '67' '100' '200'};
Periods = {'40' '20' '15' '10' '05'};
Freqs = {'25' '50' '67' '100' '200'};
N = 4;  % #of boxplot in a group
delta = linspace(-.35,.35,5); %// define offsets to distinguish plots
width = eps; %// small width to avoid overlap
legWidth = 2.3; %// make room for legend

pcell1 = [PSD40.' PSD20.' PSD15.' PSD10.' PSD5.'];
pcell2 = [(zeros(length(PSD40),1)+1).' (zeros(length(PSD20),1)+2).' ...
    (zeros(length(PSD15),1)+3).' (zeros(length(PSD10),1)+4).' (zeros(length(PSD5),1)+5).'];

gcell1 = [GSD40.' GSD20.' GSD15.' GSD10.' GSD5.'];
gcell2 = [(zeros(length(GSD40),1)+1).' (zeros(length(GSD20),1)+2).' ...
    (zeros(length(GSD15),1)+3).' (zeros(length(GSD10),1)+4).' (zeros(length(GSD5),1)+5).'];

pvel1 = [PE40 PE20 PE15 PE10 PE5];
pvel2 = [(zeros(length(PE40),1)+1).' (zeros(length(PE20),1)+2).' ...
    (zeros(length(PE15),1)+3).' (zeros(length(PE10),1)+4).' (zeros(length(PE5),1)+5).'];

gvel1 = [GE40 GE20 GE15 GE10 GE5];
gvel2 = [(zeros(length(GE40),1)+1).' (zeros(length(GE20),1)+2).' ...
    (zeros(length(GE15),1)+3).' (zeros(length(GE10),1)+4).' (zeros(length(GE5),1)+5).'];

% plot all boxplots
for ii=1:N 
    labels = Periods; 
    hold on

    if ii == 1  %% plotting Phase Cell count
        yyaxis left
        boxplot(pcell1, pcell2,'Color', DeviceColors{ii}, 'boxstyle','filled', ...
        'position',(1:numgroups)+delta(ii), 'widths',width, 'labels',labels, 'symbol', 'o', ...
        'outliersize', 1.0, 'whisker', 0,  'notch', 'off');
        medians = [nanmedian(PSD40) nanmedian(PSD20) nanmedian(PSD15) nanmedian(PSD10) nanmedian(PSD5)];
     
        pp = plot((1:numgroups)+delta(ii), medians, 'ro');
        pp.MarkerFaceColor = DeviceColors{ii};
        pp.MarkerEdgeColor = 'k';
        pp.MarkerSize = 10;

    elseif ii == 2  %% plotting Group Cell count
        yyaxis left
        boxplot(gcell1, gcell2,'Color', DeviceColors{ii}, 'boxstyle','filled', ...
        'position',(1:numgroups)+delta(ii), 'widths',width, 'labels',labels, 'symbol', 'o', ...
        'outliersize', 1.0, 'whisker', 0,  'notch', 'off');
        medians = [nanmedian(GSD40) nanmedian(GSD20) nanmedian(GSD15) nanmedian(GSD10) nanmedian(GSD5)];
    
        pp = plot((1:numgroups)+delta(ii), medians, 'ro');
        pp.MarkerFaceColor = DeviceColors{ii};
        pp.MarkerEdgeColor = 'k';
        pp.MarkerSize = 10;

    elseif ii == 3  %% plotting Phase Velocity
        yyaxis right
        boxplot(pvel1, pvel2,'Color', DeviceColors{ii}, 'boxstyle','filled', ...
        'position',(1:numgroups)+delta(ii), 'widths',width, 'labels',labels, 'symbol', 'o', ...
        'outliersize', 1.0, 'whisker', 0,  'notch', 'off');
        medians_no = [nanmedian(PE40) nanmedian(PE20) nanmedian(PE15) nanmedian(PE10) nanmedian(PE5)];
    
        pp_no = plot((1:numgroups)+delta(ii), medians_no, 'ro');
        pp_no.MarkerFaceColor = DeviceColors{ii};
        pp_no.MarkerEdgeColor = 'k';
        pp_no.MarkerSize = 10;

    elseif ii == 4  %% plotting Group Velocity
        yyaxis right
        boxplot(gvel1, gvel2,'Color', DeviceColors{ii}, 'boxstyle','filled', ...
        'position',(1:numgroups)+delta(ii), 'widths',width, 'labels',labels, 'symbol', 'o', ...
        'outliersize', 1.0, 'whisker', 0,  'notch', 'off');
        medians_no = [nanmedian(GE40) nanmedian(GE20) nanmedian(GE15) nanmedian(GE10) nanmedian(GE5)];
    
        pp_no = plot((1:numgroups)+delta(ii), medians_no, 'ro');
        pp_no.MarkerFaceColor = DeviceColors{ii};
        pp_no.MarkerEdgeColor = 'k';
        pp_no.MarkerSize = 10;
    end  
        
end

grid on;
legend(legendEntries, 'Location', 'northwest','fontsize',10);

yyaxis left
ylabel('Standard Deviation','fontsize',15);
ylim([0 2]);
% text(0.4, 1.2,'(d)','fontsize',25);
 
yyaxis right
ylabel('Percentage Error (%)','fontsize',15);
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
% figpath = '/Users/sxue3/Documents/BayMap_Figures/fig/fig5_SDErrDist.pdf';
% set(gcf,'Units','Inches');
% pos = get(gcf,'Position');
% set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
% print(gcf,figpath,'-painters', '-dpdf','-r0');