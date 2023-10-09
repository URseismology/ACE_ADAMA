%% Plot to see the observation error -- load data 
% load African country boundaries
S1 = load('/Users/sxue3/Documents/BayMap_Figures/Data/AfrCountry.mat').S1;
madgline = polyshape((S1(448).X).',(S1(448).Y).');  % load Madgascar's Big Island

% load the coast
afcoast = load('/Users/sxue3/Documents/BayMap_Figures/Data/AfricaCoast.mat');
aflon = wrapTo180(afcoast.XY(:, 1));
aflat = afcoast.XY(:,2);
coastline = polyshape(aflon, aflat);

% load the error
err_path = '/Users/sxue3/Documents/BayMap_Figures/Data/map_error/observation_errors_R8P.txt';
A = readtable(err_path,'FileType','text');
PE10= A{:,'Var8'};
senders= A{:,'Var2'};
recivers= A{:,'Var3'};
pt = A{:,'Var4'};

% load in the locations
loc_path = '/Users/sxue3/Documents/BayMap_Figures/Data/receivers.dat';
A = readtable(loc_path,'FileType','text');
lons = A{:,'Var1'};
lats = A{:,'Var2'};

%% Plot to see the observation error -- plot paths

% get the index of bad paths
badidx = [];
for i = 1:length(PE10)
    if PE10(i) > 5   % if exceeding the percentage error
        badidx = [badidx i]; 
    end
end

% get the location of the path
badsdr = senders(badidx);
badrvr = recivers(badidx);

bsdrloc = [lons(badsdr), lats(badsdr)];
brvrloc = [lons(badrvr), lats(badrvr)];

%%%%------- Start Plotting -----------%%%%
% Plot African boundary
plot(coastline, 'EdgeColor', 'k', 'LineWidth', 3, 'FaceAlpha', 0);  
hold on
% Plot Madgascar's Big Island
plot(madgline.Vertices(1:50:end, 1), madgline.Vertices(1:50:end, 2),'k', 'LineWidth', 3);

for i=1:length(badidx)
    plot([bsdrloc(i,1), brvrloc(i,1)], [bsdrloc(i,2), brvrloc(i,2)], 'Color',[0.5 0.5 0.5])
end

% plot the stations
for i = 1:length(lons)
    plot(lons(i), lats(i), 'b.')
end

title('Period 10s Err > 5%')
% set the size of the plot
x0=10;
y0=10;
fwidth=800;
fheight=800;
set(gcf,'position',[x0,y0,fwidth,fheight]);

save('./RayP8_highObsErr.mat', 'badidx', 'badsdr', 'badrvr')
%% Compute avg. velocity distribution for outliers

badpt = pt(badidx);

[baddist, ~] = distance('gc', bsdrloc(:,2),bsdrloc(:,1),brvrloc(:,2),brvrloc(:,1));
baddist = deg2km(baddist);
badvels = baddist ./badpt;

%% Plot the distribution of vel for non-outliers

% get the index of good paths
goodidx = [];
for i = 1:length(PE10)
    if PE10(i) <= 5
        goodidx = [goodidx i];
    end
end

% get the location of the path
goodsdr = senders(goodidx);
goodrvr = recivers(goodidx);

gsdrloc = [lons(goodsdr), lats(goodsdr)];
grvrloc = [lons(goodrvr), lats(goodrvr)];

goodpt = pt(goodidx);

[gooddist, ~] = distance('gc', gsdrloc(:,2),gsdrloc(:,1),grvrloc(:,2),grvrloc(:,1));
gooddist = deg2km(gooddist);
goodvels = gooddist ./goodpt;