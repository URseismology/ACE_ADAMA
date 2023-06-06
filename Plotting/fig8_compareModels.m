%%%%% Fig 8: four example locations to show inversion
%%%%% Plotting the Vs of Litho and ADAMA
%%%%% Plotting the RP & LP for ADAMA & Litho with ADAMA measurements
%%%%% Siyu Xue -- May 10. 2023
clear

%% Load data for plotting map
dsdots = load('/Users/sxue3/Documents/BayMap_Figures/Data/AfrDownSampleDots.mat');
dslat = dsdots.dslat;
dslon = dsdots.dslon;
load('/Users/sxue3/Documents/BayMap_Figures/Data/RayleighPData/plotpara.mat');
longrid = longrid(1,:);
latgrid = latgrid(:,1);

%% Load Litho model and ADAMA updates
load('/Users/sxue3/Documents/BayMap_Figures/Data/Litho1D.mat');
% load('/Users/sxue3/Documents/BayMap_Figures/Data/ADAMAinvert/CADAMA_fixed.mat');
load('/Users/sxue3/Documents/BayMap_Figures/Data/ADAMAinvert/CADAMA_fixed_new.mat')

% load the Litho velocity
load('/Users/sxue3/Documents/BayMap_Figures/Data/Litho_dwsample.mat');

% load the reference model (used to determine measurement quality)
Rayref = readtable('/Users/sxue3/Documents/BayMap_Figures/Data/SDISPR.ASC', 'FileType','fixedwidth');
Rvelref = [Rayref{78,5} Rayref{86,5} Rayref{97,5} Rayref{103,5} Rayref{108,5} ...
    Rayref{112,5} Rayref{116,5} Rayref{119,5} Rayref{121,5} Rayref{122,5} Rayref{123,5}];

%% Plotting
%%%%%%%%%%%%%%%%%%% Plot the location of the 1D model
ims = [7, 6, 28, 29];

for i = 1:4 %1:2381  % index of down sampled point
%     close all
    im = ims(i);
    if isempty(ADAMA1D{im})
        continue   % skip if no update
    end

%     % Plot the map with dots (show the location)
%     subplot(1,6,[5 6])
%     hold on
%     plot(coastline, 'EdgeColor', 'k', 'LineWidth', 2, 'FaceAlpha', 0);  
%     % Plot Madgascar's Big Island
%     plot(madgline.Vertices(1:50:end, 1), madgline.Vertices(1:50:end, 2),'k', 'LineWidth', 2);
%     % Plot African cratons
%     plot(lonCongo, latCongo, 'color', 'k', 'lineWidth', 2);
%     plot(lonKala, latKala, 'color', 'k', 'lineWidth', 2);
%     plot(lonSahara, latSahara, 'color', 'k', 'lineWidth', 2);
%     plot(lonTanza, latTanza, 'color', 'k', 'lineWidth', 2);
%     plot(lonWest, latWest, 'color', 'k', 'lineWidth', 2);
%     plot(lonZim, latZim, 'color', 'k', 'lineWidth', 2);
%     
%     % plot the dots
%     plot(dslon(im), dslat(im), '.', 'MarkerSize',15);
%     
%     % plot settings
%     minlon = -20;
%     maxlon = 55;
%     minlat = -35;
%     maxlat = 40;
%     xlim([minlon maxlon]);
%     ylim([minlat maxlat]);
%     axesm('mercator','maplatlimit',[minlat maxlat],'maplonlimit',[minlon maxlon]);
%     xticks(10*ceil(minlon/10):10:10*floor(maxlon/10));
%     yticks(10*ceil(minlat/10):10:10*floor(maxlat/10));
%     title(['point : ', num2str(im)])
%     hold off
    
    %%%%%%%%%%%%%%%%%% Plot the Litho model and ADAMA update to compare
    % convert Litho to the same format
    vSGrid = Litho1D{im}.vs;
    thckGrid = Litho1D{im}.thick;
    
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
    
    for idot = 1:length(depthGrid)
        pvs(layerst: depthGrid(idot)) = pvs(layerst: depthGrid(idot)) .* vSGrid(idot);
        layerst = 1 + depthGrid(idot);
    end
        
    % fill to the max. depth
    pvs(layerst: end) = pvs(layerst: end) .* vSGrid(end);
    
    
    % Plot Litho model with the updated model
    ADAMAvs = ADAMA1D{im}.vsv_update;
    fupdate = find(ADAMAvs(:,1) ~= 0);
    
    if i == 1
        subplot(4,4,4)
    elseif i == 2
        subplot(4,4,8)
    elseif i == 3
        subplot(4,4,12)
    else
        subplot(4,4,16)
    end

    hold on
    plot(pvs./1000, 'k', 'LineWidth',1)
    plot(ADAMAvs(fupdate(end),:)./1000, 'r', 'LineWidth',1.5)
    xlim([1 100])
    xticks([0 20 40 60 80 100])
    xticklabels({'0','10','20','30','40','50'})  % grid size is 500m
    camroll(-90)
    xlabel('Depth (km)')
    ylim([1.5 5])

    if i == 1
        ylabel('Velocity (km/s)')
    end
    
    %%%%%%%%%%%%%%%%%% Plot the measurement: Litho, ADAMA, ADAMA updates
    
    % Get the sample points in map
    outdot = [];
    afrmask = zeros(2381, 2);
    for idot = 1:length(dslat)
        xi = find(latgrid <= dslat(idot), 1, 'last');
        yi = find(longrid <= dslon(idot), 1, 'last');
    
        afrmask(idot, :) = [xi, yi];
    end
    
    % Load ADAMA velocity measurement
    type = 'P'; % 'P' or 'G'
    Periods = {'R5', 'R6', 'R8', 'R10', 'R12', 'R15', 'R20', 'R25', 'R30', 'R35', 'R40'};
    % Periods = {'L5', 'L6', 'L8', 'L10', 'L12', 'L15', 'L20', 'L25', 'L30', 'L35', 'L40'};
    
    for ip = 1:length(Periods)
        period = Periods{ip};
        
        pavgpath = ['/Users/sxue3/Documents/BayMap_Figures/Data/RayAvgMap/', period, '_', type, '.mat'];
    %     pavgpath = ['/Users/sxue3/Documents/BayMap_Figures/Data/LoveAvgMap/', period, '_', type, '_maps.mat'];
        pavg = load(pavgpath).avgmap;
        psd = load(pavgpath).sdevmap;
    
        % get the vel. and sd. at a point
        ipavg(ip) = pavg(afrmask(im, 1), afrmask(im, 2));
        ipsd(ip) = psd(afrmask(im, 1), afrmask(im, 2));
    end
    
    if i == 1
        subplot(4,4,[1 2 3])
    elseif i == 2
        subplot(4,4,[5 6 7])
    elseif i == 3
        subplot(4,4,[9 10 11])
    else
        subplot(4,4,[13 14 15])
    end

    % Plot all three velocities
    hold on
    p1 = plot(flip(LithoRay(im, :))./1000, 'b', 'LineWidth',1.5);   % Litho velocity
    p2 = plot(ADAMA1D{im}.U./1000, 'r', 'LineWidth',1.5);     % the final ADAMA update
    p3 = plot(Rvelref, 'Color',[0.502 0.502 0.502], 'LineWidth',1.5);   % the reference velocity
    p4 = errorbar(ipavg, ipsd, 'k.');    % ADAMA measurements w. uncert.
    xticks([1 2 3 4 5 6 7 8 9 10 11])
    xticklabels([5 6 8 10 12 15 20 25 30 35 40]);
    ylabel('Velocity (km/s)')
    ylim([2 6])
    hold off
    
    
    if i == 1
        legend('Litho RP', 'Final RP', 'Reference', 'Location','southeast')
    elseif i == 4
        xlabel('Period (s)')
    end


end

% set the size of the plot
x0=10;
y0=10;
fwidth=500;
fheight=600;
set(gcf,'position',[x0,y0,fwidth,fheight]);

% save the figure
figpath = '/Users/sxue3/Documents/BayMap_Figures/fig/fig8_ModelSample.pdf';
set(gcf,'Units','Inches');
pos = get(gcf,'Position');
set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
print(gcf,figpath,'-painters', '-dpdf','-r0');