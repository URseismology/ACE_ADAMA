%%%%%%% Fig suppliment: plot the boxplots of map velocity (by depth) for
%%%%%%% each crustal types (using ECM1 & Crust1.0)
%%%%%%% Siyu Xue -- AUG 20. 2023
clear

%%%% INPUT the depth of interest
mdepth = 10;    % [km] in depth

%%% NOTE: the results are not so ideal, so probably won't use this code

%% (1.1) Prepare the ADAMA Velocity map (slice)

% load ADAMA
load('/Users/sxue3/Documents/BayMap_Figures/Data/Litho1D.mat')  % load Litho 1D models
load('/Users/sxue3/Documents/BayMap_Figures/Data/ADAMAinvert/ADAMA_1D_fixed1p7.mat')  % load ADAMA inversion results
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
% Sahara
Sahara_line = polyshape(AC.Sahara(:,2), AC.Sahara(:,1));

% Congo
Congo_line = polyshape(AC.Congo(:,2), AC.Congo(:,1));

% Tanza
Tanza_line = polyshape(AC.Tanza(:,2), AC.Tanza(:,1));

% Kala
Kala_line = polyshape(AC.Kala(:,2), AC.Kala(:,1));

% West
West_line = polyshape(AC.West(:,2), AC.West(:,1));

%% Assemble Litho slice

iD = mdepth*2+1;    % depth = 500*(iD-1) [m]
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

%% (1.2) Plot the ADAMA velocity map

clim = [3.9, 4.7];  % colorbar limit

% create a colormap to show the difference
Red = [1,0,0];
White = [1,1,1];
Blue = [0.0745, 0.6235, 1];
RWBcolormap = [[ones(1,128), linspace(1,0.0745,128)]; [linspace(0, 1, 128), linspace(1,0.6235,128)]; [linspace(0, 1, 128), ones(1,128)]].';

figure(120)
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


% set the size of the plot
x0=10;
y0=10;
fwidth=700;
fheight=700;
set(gcf,'position',[x0,y0,fwidth,fheight], 'Color', 'w');
set(gca,'fontsize', 14) 

% save the figure
% figpath = '/Users/sxue3/Documents/BayMap_Figures/fig/fig11_VelGradient1.pdf';
% set(gcf,'Units','Inches');
% pos = get(gcf,'Position');
% set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
% print(gcf,figpath,'-painters', '-dpdf','-r0');


%% (2.1) Read in the Crust1.0 crustal model
C1model = readmatrix('/Users/sxue3/Documents/BayMap_Figures/Data/Crust1.0/CNtype1-1_idx.txt');

C1xgrid = -179.5:1:179.5;    % lon
C1ygrid = 89.5:-1:-89.5;    % lat
[C1xq, C1yq] = meshgrid(C1xgrid, C1ygrid);

figure(210)
pcolor(C1xq,C1yq, C1model)
shading flat
colormap(hsv)

%% (2.2) Record African crustal type at 0.1 * 0.1 degree [Crust1.0]

ADAMA_C1 = Vq_ADAMA;

for irow = 1:751
    for icol = 1:751
        if ~isnan(Vq_ADAMA(irow, icol))
            tempx = round(xq(irow, icol)-0.5)+0.5;
            tempy = round(yq(irow, icol)-0.5)+0.5;
            ADAMA_C1(irow, icol) = C1model(C1ygrid==tempy, C1xgrid==tempx);
        end
    end
end

% if ignoring the marginal crustal types
new_ADAMA_C1 = ADAMA_C1;
new_ADAMA_C1(ADAMA_C1 == 21) = NaN;
new_ADAMA_C1(ADAMA_C1 == 26) = NaN;
new_ADAMA_C1(ADAMA_C1 == 28) = NaN;
new_ADAMA_C1(ADAMA_C1 == 29) = NaN;
new_ADAMA_C1(ADAMA_C1 == 30) = NaN;
new_ADAMA_C1(ADAMA_C1 == 32) = NaN;

% create a colormap to show major African crustal types
C1cmap = ones(36,3);
C1cmap(4, :) = [0, 0.4471, 0.7412];    % early Archean
C1cmap(5, :) = [0.851, 0.3255, 0.098];    % late Archean
C1cmap(6, :) = [0.4941, 0.1843, 0.5569];    % early/mid  Proter
C1cmap(8, :) = [0.9294, 0.6941, 0.1255];    % late Proter
C1cmap(9, :) = [0.4667, 0.6745, 0.1882];    % slow late Proter
C1cmap(14, :) = [0.0588, 1, 1];   % extended crust
C1cmap(17, :) = [1, 0.0745, 0.651];   %  orogen, thick upper crust, very thin lower crust
C1cmap(23, :) = [0.502, 0.502, 0.502];   % Rift
C1cmap(25, :) = [0, 0, 1];   % fast Phanerozoic (E. Australia, S. Africa, N. Siberia)

figure(220)
pcolor(xq,yq, new_ADAMA_C1)
shading flat
colormap(C1cmap)
caxis([1, 36])

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
title('Crust1.0 Model')
set(gca,'fontsize', 14, 'Color', 'w') 



%% manually create a legend for crustal types
figure(221)
bar(0,1,'FaceColor', C1cmap(4, :))
hold on
bar(0,1,'FaceColor', C1cmap(5, :))
bar(0,1,'FaceColor', C1cmap(6, :))
bar(0,1,'FaceColor', C1cmap(8, :))
bar(0,1,'FaceColor', C1cmap(9, :))
bar(0,1,'FaceColor', C1cmap(14, :))
bar(0,1,'FaceColor', C1cmap(17, :))
bar(0,1,'FaceColor', C1cmap(23, :))
bar(0,1,'FaceColor', C1cmap(25, :))

legend({'early Archean','late Archean','early/mid  Proter','late  Proter','slow late Proter',...
    'extended crust','orogen','Rift','fast Phanerozoic'})

%% (3.1) Read in the ECM1 crustal model
ECM1 = readtable('/Users/sxue3/Documents/BayMap_Figures/Data/ECM1-CrustalTypes/CrustMapLayers_idx.txt');

EClon = ECM1{:,'Lon'};
EClat = ECM1{:,'Lat'};
ECtype = ECM1{:,'Type'};

[ECxg,ECyg] = meshgrid(unique(EClon), unique(EClat));
F = scatteredInterpolant(EClon, EClat, ECtype);
ECmodel = F(ECxg,ECyg);

figure(310)
pcolor(ECxg,ECyg, ECmodel)
shading flat
colormap(turbo)
colorbar

%% (3.2) Record African crustal type at 0.1 * 0.1 degree [ECM1]

ADAMA_EC = Vq_ADAMA;

unilat = unique(EClat);
unilon = unique(EClon);

for irow = 1:751
    for icol = 1:751
        if ~isnan(Vq_ADAMA(irow, icol))
            tempx = round(xq(irow, icol)-0.5)+0.5;
            tempy = round(yq(irow, icol)-0.5)+0.5;
            ADAMA_EC(irow, icol) = ECmodel(unilat==tempy, unilon==tempx);
        end
    end
end

% if ignoring the marginal crustal types
new_ADAMA_EC = ADAMA_EC;
new_ADAMA_EC(ADAMA_EC == 1) = NaN;
new_ADAMA_EC(ADAMA_EC == 2) = NaN;
new_ADAMA_EC(ADAMA_EC == 3) = NaN;
new_ADAMA_EC(ADAMA_EC == 4) = NaN;

% create a colormap to show major African crustal types
ECcmap = ones(9,3);
ECcmap(5, :) = [0.851, 0.3255, 0.098];    % Extended Crust
ECcmap(6, :) = [0.9294, 0.6941, 0.1255];    % Shield
ECcmap(7, :) = [0.4941, 0.1843, 0.5569];    % Orogen
ECcmap(8, :) = [0, 0.4471, 0.7412];    % Platform
ECcmap(9, :) = [0.4667, 0.6745, 0.1882];    % Basin

figure(320)
pcolor(xq,yq, new_ADAMA_EC)
shading flat
colormap(ECcmap)
caxis([1 9])

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
title('ECM1 Model')
set(gca,'fontsize', 14, 'Color', 'w') 


%% manually create a legend for crustal types
figure(321)
bar(0,1,'FaceColor', ECcmap(5, :))
hold on
bar(0,1,'FaceColor', ECcmap(6, :))
bar(0,1,'FaceColor', ECcmap(7, :))
bar(0,1,'FaceColor', ECcmap(8, :))
bar(0,1,'FaceColor', ECcmap(9, :))
legend({'Extended Crust','Shield','Orogen','Platform','Basin'})

%% (4.1) Plot the crustal velocity boxplot using Crust1.0

% plot the Crust1.0 boxplot (all types)
% boxplot(reshape(Vq_ADAMA/1000,[], 1), reshape(ADAMA_C1,[], 1), 'Labels',...
%     {'early Archean','late Archean','early/mid  Proter','late  Proter','slow late Proter',...
%     'extended crust','orogen','Margin-continent/shield  transition','Rift','fast Phanerozoic',...
%     'normal oceanic','oceanic plateaus','continental shelf','continental slope','thinned cont. crust'})

figure(410)
% plot the Crust1.0 boxplot (only major types: 4,5,6,8,9,14,17,23,25)
boxplot(reshape(Vq_ADAMA/1000,[], 1), reshape(new_ADAMA_C1,[], 1), 'Labels',...
    {'early Archean','late Archean','early/mid  Proter','late  Proter','slow late Proter',...
    'extended crust','orogen','Rift','fast Phanerozoic'})
title(['Crust1.0 - ', num2str(mdepth), 'km'])

%% (4.2) Plot the crustal velocity boxplot using ECM1

figure(420)
% plot only the major types: 5,6,7,8,9
boxplot(reshape(Vq_ADAMA/1000,[], 1), reshape(new_ADAMA_EC,[], 1), 'Labels',...
    {'Extended Crust','Shield','Orogen','Platform','Basin'})
title(['ECM1 - ', num2str(mdepth), 'km'])

%% (5.1) Plot the crustal velocity boxplot for individual craton using Crust1.0

C1_sahara = new_ADAMA_C1;
C1_west = new_ADAMA_C1;
C1_kala = new_ADAMA_C1;
C1_congo = new_ADAMA_C1;
C1_tanza = new_ADAMA_C1;

for irow = 1:751
    for icol = 1:751
        if ~isnan(new_ADAMA_C1(irow, icol))
            if ~isinterior(Sahara_line, [xq(irow, icol), yq(irow, icol)])
                C1_sahara(irow, icol) = NaN;
            end
            if ~isinterior(West_line, [xq(irow, icol), yq(irow, icol)])
                C1_west(irow, icol) = NaN;
            end
            if ~isinterior(Kala_line, [xq(irow, icol), yq(irow, icol)])
                C1_kala(irow, icol) = NaN;
            end
            if ~isinterior(Congo_line, [xq(irow, icol), yq(irow, icol)])
                C1_congo(irow, icol) = NaN;
            end
            if ~isinterior(Tanza_line, [xq(irow, icol), yq(irow, icol)])
                C1_tanza(irow, icol) = NaN;
            end
        end
    end
end


figure(510)

subplot(2,3,1)  % Sahara
boxplot(reshape(Vq_ADAMA(~isnan(C1_sahara))/1000,[], 1), reshape(C1_sahara(~isnan(C1_sahara)),[], 1), 'Labels',...
    {'early Archean','late  Proter','slow late Proter','extended crust'})
title(['Crust1.0 - Sahara - ', num2str(mdepth), 'km'])

subplot(2,3,2)  % West
boxplot(reshape(Vq_ADAMA(~isnan(C1_west))/1000,[], 1), reshape(C1_west(~isnan(C1_west)),[], 1), 'Labels',...
    {'early Archean','late Archean','early/mid  Proter','late  Proter','orogen'})
title(['Crust1.0 - West - ', num2str(mdepth), 'km'])

subplot(2,3,3)  % Kala
boxplot(reshape(Vq_ADAMA(~isnan(C1_kala))/1000,[], 1), reshape(C1_kala(~isnan(C1_kala)),[], 1), 'Labels',...
    {'early Archean','early/mid  Proter','late  Proter','fast Phanerozoic'})
title(['Crust1.0 - Kala - ', num2str(mdepth), 'km'])

subplot(2,3,4)  % Congo
boxplot(reshape(Vq_ADAMA(~isnan(C1_congo))/1000,[], 1), reshape(C1_congo(~isnan(C1_congo)),[], 1), 'Labels',...
    {'early Archean','late Archean','early/mid  Proter','late  Proter','Rift'})
title(['Crust1.0 - Congo - ', num2str(mdepth), 'km'])

subplot(2,3,5)  % Tanza
boxplot(reshape(Vq_ADAMA(~isnan(C1_tanza))/1000,[], 1), reshape(C1_tanza(~isnan(C1_tanza)),[], 1), 'Labels',...
    {'early Archean','late Archean','early/mid  Proter','late  Proter','Rift'})
title(['Crust1.0 - Tanza - ', num2str(mdepth), 'km'])


%% (6.1) Plot the crustal velocity boxplot for individual craton using ECM1

EC_sahara = new_ADAMA_EC;
EC_west = new_ADAMA_EC;
EC_kala = new_ADAMA_EC;
EC_congo = new_ADAMA_EC;
EC_tanza = new_ADAMA_EC;

for irow = 1:751
    for icol = 1:751
        if ~isnan(new_ADAMA_EC(irow, icol))
            if ~isinterior(Sahara_line, [xq(irow, icol), yq(irow, icol)])
                EC_sahara(irow, icol) = NaN;
            end
            if ~isinterior(West_line, [xq(irow, icol), yq(irow, icol)])
                EC_west(irow, icol) = NaN;
            end
            if ~isinterior(Kala_line, [xq(irow, icol), yq(irow, icol)])
                EC_kala(irow, icol) = NaN;
            end
            if ~isinterior(Congo_line, [xq(irow, icol), yq(irow, icol)])
                EC_congo(irow, icol) = NaN;
            end
            if ~isinterior(Tanza_line, [xq(irow, icol), yq(irow, icol)])
                EC_tanza(irow, icol) = NaN;
            end
        end
    end
end

% plot only the major types: 5,6,7,8,9
% {'Extended Crust','Shield','Orogen','Platform','Basin'}

figure(610)

subplot(2,3,1)  % Sahara
boxplot(reshape(Vq_ADAMA(~isnan(EC_sahara))/1000,[], 1), reshape(EC_sahara(~isnan(EC_sahara)),[], 1)...
    ,'Labels',{'Extended Crust','Shield','Platform'})
title(['ECM1 - Sahara - ', num2str(mdepth), 'km'])

subplot(2,3,2)  % West
boxplot(reshape(Vq_ADAMA(~isnan(EC_west))/1000,[], 1), reshape(EC_west(~isnan(EC_west)),[], 1)...
    ,'Labels',{'Extended Crust','Shield','Platform','Basin'})
title(['ECM1 - West - ', num2str(mdepth), 'km'])

subplot(2,3,3)  % Kala
boxplot(reshape(Vq_ADAMA(~isnan(EC_kala))/1000,[], 1), reshape(EC_kala(~isnan(EC_kala)),[], 1)...
    ,'Labels',{'Extended Crust','Shield','Platform','Basin'})
title(['ECM1 - Kala - ', num2str(mdepth), 'km'])

subplot(2,3,4)  % Congo
boxplot(reshape(Vq_ADAMA(~isnan(EC_congo))/1000,[], 1), reshape(EC_congo(~isnan(EC_congo)),[], 1)...
    ,'Labels',{'Extended Crust','Shield','Platform','Basin'})
title(['ECM1 - Congo - ', num2str(mdepth), 'km'])

subplot(2,3,5)  % Tanza
boxplot(reshape(Vq_ADAMA(~isnan(EC_tanza))/1000,[], 1), reshape(EC_tanza(~isnan(EC_tanza)),[], 1)...
    ,'Labels',{'Extended Crust','Shield'})
title(['ECM1 - Tanza - ', num2str(mdepth), 'km'])
