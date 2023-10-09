%%%% Plot the boxplot to show the distribution of cell and velocity in Africa
%%%% (not included in the paper)
%%%%% Siyu Xue -- Feb 10. 2023

addpath('/Users/sxue3/Documents/BayMap_Figures/Utils')
%% Read in the data
PC40 = load('/Users/sxue3/Documents/BayMap_Figures/Data/CellVelDistribution/R40_P.mat').cellinafr;
PV40 = load('/Users/sxue3/Documents/BayMap_Figures/Data/CellVelDistribution/R40_P.mat').velinafr;
GC40 = load('/Users/sxue3/Documents/BayMap_Figures/Data/CellVelDistribution/R40_G.mat').cellinafr;
GV40 = load('/Users/sxue3/Documents/BayMap_Figures/Data/CellVelDistribution/R40_G.mat').velinafr;
PC20 = load('/Users/sxue3/Documents/BayMap_Figures/Data/CellVelDistribution/R20_P.mat').cellinafr;
PV20 = load('/Users/sxue3/Documents/BayMap_Figures/Data/CellVelDistribution/R20_P.mat').velinafr;
GC20 = load('/Users/sxue3/Documents/BayMap_Figures/Data/CellVelDistribution/R20_G.mat').cellinafr;
GV20 = load('/Users/sxue3/Documents/BayMap_Figures/Data/CellVelDistribution/R20_G.mat').velinafr;
PC15 = load('/Users/sxue3/Documents/BayMap_Figures/Data/CellVelDistribution/R15_P.mat').cellinafr;
PV15 = load('/Users/sxue3/Documents/BayMap_Figures/Data/CellVelDistribution/R15_P.mat').velinafr;
GC15 = load('/Users/sxue3/Documents/BayMap_Figures/Data/CellVelDistribution/R15_G.mat').cellinafr;
GV15 = load('/Users/sxue3/Documents/BayMap_Figures/Data/CellVelDistribution/R15_G.mat').velinafr;
PC10 = load('/Users/sxue3/Documents/BayMap_Figures/Data/CellVelDistribution/R10_P.mat').cellinafr;
PV10 = load('/Users/sxue3/Documents/BayMap_Figures/Data/CellVelDistribution/R10_P.mat').velinafr;
GC10 = load('/Users/sxue3/Documents/BayMap_Figures/Data/CellVelDistribution/R10_G.mat').cellinafr;
GV10 = load('/Users/sxue3/Documents/BayMap_Figures/Data/CellVelDistribution/R10_G.mat').velinafr;
PC5 = load('/Users/sxue3/Documents/BayMap_Figures/Data/CellVelDistribution/R5_P.mat').cellinafr;
PV5 = load('/Users/sxue3/Documents/BayMap_Figures/Data/CellVelDistribution/R5_P.mat').velinafr;
GC5 = load('/Users/sxue3/Documents/BayMap_Figures/Data/CellVelDistribution/R5_G.mat').cellinafr;
GV5 = load('/Users/sxue3/Documents/BayMap_Figures/Data/CellVelDistribution/R5_G.mat').velinafr;

%% initialize some parameters
numgroups = 5;
DeviceColors = {[0.59 0.89 1],[0 0.45 0.74],[1 0.82 0.46], [1 0.41 0.01]};
legendEntries = {'P-Cell','G-Cell', 'P-Vel','G-Vel'};
% Periods = {'40' '20' '15' '10' '05'};
% Freqs = {'25' '50' '67' '100' '200'};
Periods = {'40' '20' '15' '10' '05'};
Freqs = {'25' '50' '67' '100' '200'};
N = 4;  % #of boxplot in a group
delta = linspace(-.35,.35,5); %// define offsets to distinguish plots
width = eps; %// small width to avoid overlap
legWidth = 2.3; %// make room for legend

pcell1 = [PC40.' PC20.' PC15.' PC10.' PC5.'];
pcell2 = [(zeros(length(PC40),1)+1).' (zeros(length(PC20),1)+2).' ...
    (zeros(length(PC15),1)+3).' (zeros(length(PC10),1)+4).' (zeros(length(PC5),1)+5).'];

gcell1 = [GC40.' GC20.' GC15.' GC10.' GC5.'];
gcell2 = [(zeros(length(GC40),1)+1).' (zeros(length(GC20),1)+2).' ...
    (zeros(length(GC15),1)+3).' (zeros(length(GC10),1)+4).' (zeros(length(GC5),1)+5).'];

pvel1 = [PV40 PV20 PV15 PV10 PV5];
pvel2 = [(zeros(length(PV40),1)+1).' (zeros(length(PV20),1)+2).' ...
    (zeros(length(PV15),1)+3).' (zeros(length(PV10),1)+4).' (zeros(length(PV5),1)+5).'];

gvel1 = [GV40 GV20 GV15 GV10 GV5];
gvel2 = [(zeros(length(GV40),1)+1).' (zeros(length(GV20),1)+2).' ...
    (zeros(length(GV15),1)+3).' (zeros(length(GV10),1)+4).' (zeros(length(GV5),1)+5).'];

% plot all boxplots
for ii=1:N 
    labels = Periods; 
    hold on

    if ii == 1  %% plotting Phase Cell count
        yyaxis left
        boxplot(pcell1, pcell2,'Color', DeviceColors{ii}, 'boxstyle','filled', ...
        'position',(1:numgroups)+delta(ii), 'widths',width, 'labels',labels, 'symbol', 'o', ...
        'outliersize', 1.0, 'whisker', 0,  'notch', 'off');
        medians = [nanmedian(PC40) nanmedian(PC20) nanmedian(PC15) nanmedian(PC10) nanmedian(PC5)];
     
        pp = plot((1:numgroups)+delta(ii), medians, 'ro');
        pp.MarkerFaceColor = DeviceColors{ii};
        pp.MarkerEdgeColor = 'k';
        pp.MarkerSize = 10;

    elseif ii == 2  %% plotting Group Cell count
        yyaxis left
        boxplot(gcell1, gcell2,'Color', DeviceColors{ii}, 'boxstyle','filled', ...
        'position',(1:numgroups)+delta(ii), 'widths',width, 'labels',labels, 'symbol', 'o', ...
        'outliersize', 1.0, 'whisker', 0,  'notch', 'off');
        medians = [nanmedian(GC40) nanmedian(GC20) nanmedian(GC15) nanmedian(GC10) nanmedian(GC5)];
    
        pp = plot((1:numgroups)+delta(ii), medians, 'ro');
        pp.MarkerFaceColor = DeviceColors{ii};
        pp.MarkerEdgeColor = 'k';
        pp.MarkerSize = 10;

    elseif ii == 3  %% plotting Phase Velocity
        yyaxis right
        boxplot(pvel1, pvel2,'Color', DeviceColors{ii}, 'boxstyle','filled', ...
        'position',(1:numgroups)+delta(ii), 'widths',width, 'labels',labels, 'symbol', 'o', ...
        'outliersize', 1.0, 'whisker', 0,  'notch', 'off');
        medians_no = [nanmedian(PV40) nanmedian(PV20) nanmedian(PV15) nanmedian(PV10) nanmedian(PV5)];
    
        pp_no = plot((1:numgroups)+delta(ii), medians_no, 'ro');
        pp_no.MarkerFaceColor = DeviceColors{ii};
        pp_no.MarkerEdgeColor = 'k';
        pp_no.MarkerSize = 10;

    elseif ii == 4  %% plotting Group Velocity
        yyaxis right
        boxplot(gvel1, gvel2,'Color', DeviceColors{ii}, 'boxstyle','filled', ...
        'position',(1:numgroups)+delta(ii), 'widths',width, 'labels',labels, 'symbol', 'o', ...
        'outliersize', 1.0, 'whisker', 0,  'notch', 'off');
        medians_no = [nanmedian(GV40) nanmedian(GV20) nanmedian(GV15) nanmedian(GV10) nanmedian(GV5)];
    
        pp_no = plot((1:numgroups)+delta(ii), medians_no, 'ro');
        pp_no.MarkerFaceColor = DeviceColors{ii};
        pp_no.MarkerEdgeColor = 'k';
        pp_no.MarkerSize = 10;
    end  
        
end

grid on;
legend(legendEntries, 'Location', 'northwest','fontsize',10);

yyaxis left
ylabel('Cell Count','fontsize',15);


% ylim([0 7]);
% text(0.4, 1.2,'(d)','fontsize',25);
 
yyaxis right
ylabel('Velocity (km/s)','fontsize',15);
% ylim([-3 3]);

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
% figpath = '/Users/sxue3/Documents/BayMap_Figures/fig/fig3_RayPGDist.pdf';
% set(gcf,'Units','Inches');
% pos = get(gcf,'Position');
% set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
% print(gcf,figpath,'-painters', '-dpdf','-r0');