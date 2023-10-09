%%%%% Fig 9: transec maps of Lttho, ADAMA, Litho-ADAMA difference, and uncertanities
%%%%% also the topography along the transec
%%%%% Siyu Xue -- May 10, 2023

addpath('/Users/sxue3/Documents/ADAMA_Figures/m_map')

%% (1) Read in data
load('/Users/sxue3/Documents/BayMap_Figures/Data/Litho1D.mat')  % load Litho 1D models
load('/Users/sxue3/Documents/BayMap_Figures/Data/ADAMAinvert/CADAMA_fixed.mat')  % load ADAMA inversion results
load('/Users/sxue3/Documents/BayMap_Figures/Data/AfrDownSampleDots.mat')  % load data locations
% read in the moho data file
A = readmatrix('/Users/sxue3/Documents/BayMap_Figures/Data/Africa_Moho_density_contrast.txt');

%% Select points based on upper and lower bounds
tslat = [];
tslon = [];
tsi = [];

%%%%%%%%%%%%%%%%%%%%%%%%%% Case (1) %%%%%%%%%%%%%%%%%%%%%%%%%%%
slope = 24/25;
b1 = -21;       % for upper bound
b2 = -22.5;     % for lower bound

for ip = 1:length(dslat)
    if dslat(ip) <= (dslon(ip)*slope + b1) && dslat(ip) >= (dslon(ip)*slope + b2)
        tslon = [tslon, dslon(ip)];
        tslat = [tslat, dslat(ip)];
        tsi = [tsi, ip];
    end
end

% manually add/delete points
tslon(5) = [];
tslat(5) = [];
tsi(5) = [];

% Use this online calculator: https://keisan.casio.com/exec/system/1224587128
totalDist = 4015;       % total distance between the two end points

% Sort all points based on increasing lon
[tslon, I] = sort(tslon);
tslat = tslat(I);
tsi = tsi(I);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %%%%%%%%%%%%%%%%%%%%%%%%%%% Case (2) %%%%%%%%%%%%%%%%%%%%%%%%%%%
% slope = 4/10;
% b1 = -30;       % for upper bound
% b2 = -31.2;     % for lower bound
% 
% for ip = 1:length(dslat)
%     if dslat(ip) <= (dslon(ip)*slope + b1) && dslat(ip) >= (dslon(ip)*slope + b2)
%         tslon = [tslon, dslon(ip)];
%         tslat = [tslat, dslat(ip)];
%         tsi = [tsi, ip];
%     end
% end
% 
% % Use this online calculator: https://keisan.casio.com/exec/system/1224587128
% totalDist =2777;       % total distance between the two end points
% 
% % Sort all points based on increasing lon
% [tslon, I] = sort(tslon);
% tslat = tslat(I);
% tsi = tsi(I);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %%%%%%%%%%%%%%%%%%%%%%%%%%% Case (3) %%%%%%%%%%%%%%%%%%%%%%%%%%%
% slope = -1/2;
% b1 = 27;       % for upper bound
% b2 = 25.8;     % for lower bound
% 
% for ip = 1:length(dslat)
%     if dslat(ip) <= (dslon(ip)*slope + b1) && dslat(ip) >= (dslon(ip)*slope + b2)
%         tslon = [tslon, dslon(ip)];
%         tslat = [tslat, dslat(ip)];
%         tsi = [tsi, ip];
%     end
% end
% 
% % Use this online calculator: https://keisan.casio.com/exec/system/1224587128
% totalDist = 6586;       % total distance between the two end points
% 
% % Sort all points based on increasing lon
% [tslon, I] = sort(tslon);
% tslat = tslat(I);
% tsi = tsi(I);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %%%%%%%%%%%%%%%%%%%%%%%%%%% Case (4) %%%%%%%%%%%%%%%%%%%%%%%%%%%
% slope = -4;
% b1 = 86.2;     % for upper bound
% b2 = 82.2;     % for lower bound
% 
% for ip = 1:length(dslat)
%     if dslat(ip) <= (dslon(ip)*slope + b1) && dslat(ip) >= (dslon(ip)*slope + b2)
%         tslon = [tslon, dslon(ip)];
%         tslat = [tslat, dslat(ip)];
%         tsi = [tsi, ip];
%     end
% end
% 
% % manually add/delete points
% addpoints = [2262, 1070, 885, 1640];
% tsi = [tsi, addpoints];
% tslat = [tslat, dslat(addpoints).'];
% tslon = [tslon, dslon(addpoints).'];
% 
% % Use this online calculator: https://keisan.casio.com/exec/system/1224587128
% totalDist = 7292;       % total distance between the two end points
% 
% % Sort all points based on increasing lon
% [tslat, I] = sort(tslat);
% tslon = tslon(I);
% tsi = tsi(I);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% get the moho lat and lon
moho_lat = round(tslat);
moho_lon = round(tslon);

% % Plot to check points selected
% subplot(2,4,[1,2,5,6])
% plot(dslon, dslat, 'b.')
% hold on
% plot(tslon, tslat, 'r.')
% xs = -20:0.1:55;
% ys1 = slope * xs + b1;
% ys2 = slope * xs + b2;
% 
% % plot cratons
% plot(lonCongo, latCongo, 'color', 'k', 'lineWidth', 2);
% plot(lonKala, latKala, 'color', 'k', 'lineWidth', 2);
% plot(lonSahara, latSahara, 'color', 'k', 'lineWidth', 2);
% plot(lonTanza, latTanza, 'color', 'k', 'lineWidth', 2);
% plot(lonWest, latWest, 'color', 'k', 'lineWidth', 2);
% plot(lonZim, latZim, 'color', 'k', 'lineWidth', 2);

% % plot the sample selection range
% plot(xs, ys1, 'g')
% plot(xs, ys2, 'g')
% % plot the moho selection
% plot(moho_lon, moho_lat, 'k')
% xlabel('Longitude')
% ylabel('Latitude')
% minlon = -20;
% maxlon = 55;
% minlat = -35;
% maxlat = 40;
% xlim([minlon maxlon]);
% ylim([minlat maxlat]);
% axesm('mercator','maplatlimit',[minlat maxlat],'maplonlimit',[minlon maxlon]);
% xticks(10*ceil(minlon/10):10:10*floor(maxlon/10));
% yticks(10*ceil(minlat/10):10:10*floor(maxlat/10));

%% Assemble 1-D models into 2-D

% assemble the updated model
vel2D = zeros(600, length(tsi));
icol = 1;
for ip = tsi
    if ~isempty(ADAMA1D{ip})  % use ADAMA if there is an update
        ADAMAvs = ADAMA1D{ip}.vsv_update;
        fupdate = find(ADAMAvs(:,1) ~= 0);
        vel2D(:,icol) = ADAMAvs(fupdate(end),:);
    else
        Lithovs = Litho1D{ip}.pvs;
        vel2D(:,icol) = Lithovs(1:600);
    end
    icol = icol + 1;
end


%% Plot the updated model
xs = zeros(120, length(tsi));
ys = zeros(120, length(tsi));

for i = 1:120
    xs(i, :) = linspace(1,totalDist,length(tsi));
    ys(i, :) = repelem(i*0.5,length(tsi));
end

ygrid = 0.5:0.5:60;     % grid depth = 500m
xgrid = 1:1:totalDist;
[xq, yq] = meshgrid(xgrid, ygrid);

vel2D_top = vel2D(1:120, :);
Vq = interp2(xs,ys,vel2D_top,xq,yq,'linear');

% get the homo depth (with sediments) of the cross-section
moho_depth = zeros(size(moho_lon));
for im = 1:length(moho_lon)
    irow = (A(:,1) == moho_lon(im)) .* (A(:,2) == moho_lat(im));
    moho_depth(im) = A(logical(irow), 4);
end

% create a colormap to show the velocity

color1 = [1 0 0];
color13 = [1 0.15 0];
color14 = [1 0.3 0];
color15 = [1 0.45 0];
color16 = [1 0.6549 0];
color2 = [1 0.85 0];
color25 = [0.9333 1 0];
color3 = [0.8 1 0];
color4 = [0.1333 1 0];
color5 = [0 0.8667 1];
color55 = [0 0.6 1];
color6 = [0 0.4667 1];
color7 = [0.333 0 1];

allcolors = [color1; color13; color14; color15; color16; color2; color25; color3; color4; color5; color55; color6; color7];

for ic = 1:length(allcolors)-1

    istart = (ic-1)*10+1;
    iend = ic*10;

    temp = [linspace(allcolors(ic, 1), allcolors(ic+1, 1), 10); ...
        linspace(allcolors(ic, 2), allcolors(ic+1, 2), 10); linspace(allcolors(ic, 3), allcolors(ic+1, 3), 10)];
    vel_colmap(istart:iend, 1:3) = temp';
end

% plot the updated model - ADAMA
ax1 = subplot(2,4,[1,2]);
% load('./Utils/tolu_stdmap2.mat')    % load the colormap
imagesc(Vq/1000)
yticklabels({'10','20','30','40','50','60'})
hold on
plot(linspace(1,totalDist, length(tslon)), moho_depth*2, 'k', 'LineWidth', 2)  % convert moho from (km) to [500m]
% ch1 = colorbar;
% ch1.Title.String = 'v_s (km/s)';
colormap(ax1, vel_colmap)
caxis([2.5 4.8])
% xlabel('Distance [m]')
ylabel('Depth (km)')
title('ADAMA Velocity Model')

%% assemble the Litho1.0 model
vel2DLT = zeros(600, length(tsi));
icol = 1;
for ip = tsi
    Lithovs = Litho1D{ip}.pvs;
    vel2DLT(:,icol) = Lithovs(1:600);
    icol = icol + 1;
end

vel2DLT_top = vel2DLT(1:120, :);
VqLT = interp2(xs,ys,vel2DLT_top,xq,yq, 'linear');

% Plot the Litho model
ax2 = subplot(2,4,[3,4]);
imagesc(VqLT/1000)
caxis([2.5 4.8])
yticklabels({'10','20','30','40','50','60'})
ch2 = colorbar;
ch2.Title.String = 'v_s (km/s)';
colormap(ax2, vel_colmap)
xlabel('Distance (km)')
ylabel('Depth (km)')
title('Litho1.0 Velocity Model')

%% Plot the updated value

ax3 = subplot(2,4,[5,6]);
Vq_update = Vq - VqLT;
imagesc(Vq_update/1000)
yticklabels({'10','20','30','40','50','60'})
% ch3 = colorbar;
% ch3.Title.String = 'v_s (km/s)';
colormap(ax3, hot)
caxis([0 1.2])
xlabel('Distance (km)')
ylabel('Depth (km)')
title('ADAMA - Litho1.0')


% % Plot the Litho1.0 model
% figure(101)
% imagesc(VqLT)
% yticklabels({'10','20','30','40','50','60'})
% ch = colorbar;
% ch.Title.String = 'v_s [m/s]';
% colormap(parula)
% xlabel('Distance (km)')
% ylabel('Depth (km)')
% title('Litho1.0 Model')

%% Plot the ADAMA Uncertainty

Vq_uncer = NaN(120, length(tsi));

% loop through all points
for idot = 1:length(tsi)
    simu_path = ['/Users/sxue3/Documents/BayMap_Figures/Data/TransecUncer/' num2str(tsi(idot)), '.mat'];

    if exist(simu_path, "file")
    % read in the simulations if ADAMA updated
        load(simu_path)

    ADAMAidot = [];
    % save all simulations to a matrix
    for im = 1:100
        if ~isempty(ADAMAid{im})
            ADAMAvs = ADAMAid{im}.vsv_update;
            fupdate = find(ADAMAvs(:,1) ~= 0);
            ADAMAidot = [ADAMAidot; ADAMAvs(fupdate(end),1:120)];
        end 
    end

    % Compute the SD at each 500m depth
    Vq_uncer(:, idot) = std(ADAMAidot, 0, 1);
    else
        continue
    end
end

ax4 = subplot(2,4,[7,8]);
imagesc(Vq_uncer/1000)
yticklabels({'10','20','30','40','50','60'})
ch4 = colorbar;
ch4.Title.String = 'v_s (km/s)';
colormap(ax4, hot)
xlabel('Distance (km)')
ylabel('Depth (km)')
title('ADAMA Uncertainty')
caxis([0 1.2])


% set the size of the plot
x0=10;
y0=10;
fwidth=1000;
fheight=500;
set(gcf,'position',[x0,y0,fwidth,fheight]);

% save the figure
fig = gcf;
% saveFig('fig9_TransecMaps.pdf', '/Users/sxue3/Documents/BayMap_Figures/fig/', 1, fig);

%% Plot topography of the transec

% get the topo along the transec
[ELEV,LONG,LAT]=m_etopo2([min(tslon) max(tslon) min(tslat) max(tslat)]);
topo_trasec = diag(flip(ELEV));

% read in the inversion velocity quality
load('./Data/RayVelQly.mat')    
qly_trasec = dp_record(tsi);
qly_xx = linspace(1, length(topo_trasec), length(tsi));


% plot the topo
figure(22)
plot(smooth(topo_trasec, 5)/1000, 'k', 'LineWidth',1.5)
hold on

for idot = 1:length(tsi)
    if qly_trasec(idot) == 1
        plot(qly_xx(idot), 0, '.', 'MarkerSize', 20, 'Color',[0.1059 0.9804 0.1059])
    elseif qly_trasec(idot) == 2
        plot(qly_xx(idot), 0, '.', 'MarkerSize', 20, 'Color',[0.3020 0.7451 0.9333])
    elseif qly_trasec(idot) == 3
        plot(qly_xx(idot), 0, '.', 'MarkerSize', 20, 'Color',[1 0 0])
    else
        plot(qly_xx(idot), 0, '.', 'MarkerSize', 20, 'Color',[0.6353 0.0784 0.1843])
    end
end

ylabel('H (km)')
xlim([1 1535])

% set the size of the plot
x0=10;
y0=10;
fwidth=450;
fheight=60;
set(gcf,'position',[x0,y0,fwidth,fheight]);

% save the figure
fig = gcf;
% saveFig('fig9_TransecTopo.pdf', '/Users/sxue3/Documents/BayMap_Figures/fig/', 1, fig);