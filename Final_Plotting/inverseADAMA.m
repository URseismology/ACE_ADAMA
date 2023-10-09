%%%%% inverse ADAMA data using Raylee inversion code at downsampled locations 
addpath(genpath('/Users/sxue3/Documents/Hawaii/Raylee_n_Lovee/raylee_codes-master'))
addpath('/Users/sxue3/Documents/BayMap_Figures/Utils')

load('/Users/sxue3/Documents/BayMap_Figures/Data/RayleighPData/plotpara.mat');
longrid = longrid(1,:);
latgrid = latgrid(:,1);
% % load the coast
% afcoast = load('/Users/sxue3/Documents/BayMap_Figures/Data/AfricaCoast.mat');
% aflon = wrapTo180(afcoast.XY(:, 1));
% aflat = afcoast.XY(:,2);
% coastline = polyshape(aflon, aflat);
% % load Madgascar's Big Island
% S1 = load('/Users/sxue3/Documents/BayMap_Figures/Data/AfrCountry.mat').S1;
% madgline = polyshape((S1(448).X).',(S1(448).Y).');  

% Get the downsampled velocity
dsdots = load('/Users/sxue3/Documents/BayMap_Figures/Data/AfrDownSampleDots.mat');
dslat = dsdots.dslat;
dslon = dsdots.dslon;

% Get the sample points in map
afrmask = zeros(2380, 2);
for i = 1:length(dslat)
    xi = find(latgrid <= dslat(i), 1, 'last');
    yi = find(longrid <= dslon(i), 1, 'last');

    afrmask(i, :) = [xi, yi];
end

type = 'P'; % 'P' or 'G'
Periods = {'R5', 'R6', 'R8', 'R10', 'R12', 'R15', 'R20', 'R25', 'R30', 'R35', 'R40'};
% Periods = {'L5', 'L6', 'L8', 'L10', 'L12', 'L15', 'L20', 'L25', 'L30', 'L35', 'L40'};

for ip = 1:length(Periods)
    period = Periods{ip};
    
    pavgpath = ['/Users/sxue3/Documents/BayMap_Figures/Data/RayAvgMap/', period, '_', type, '_maps.mat'];
%     pavgpath = ['/Users/sxue3/Documents/BayMap_Figures/Data/LoveAvgMap/', period, '_', type, '_maps.mat'];
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

% change the scale of ADAMA uncertanities
% psds = psds .* 0.5;

%% Inverse ADAMA using Litho as the starting model
% (1) %%%%%%%%%%%%%%%%%%%%%%%%%%%%% Load data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('/Users/sxue3/Documents/BayMap_Figures/Data/Litho1D.mat')


ADAMA1D = cell(length(dslon), 1);
notinverse = [];
 
ims = [7, 6, 28, 29];     % Siyu: for testing
for im = ims %1:length(dslon)
    disp(['model: ', num2str(im)])
    vSGrid = Litho1D{im}.vs;
    vPGrid = Litho1D{im}.vp;
    rhoGrid = Litho1D{im}.rho;
    thckGrid = Litho1D{im}.thick;
    
    % in case we get fluid layer, delete it
    if vSGrid(1) == 0
        vPGrid(1) = [];
        rhoGrid(1) = [];
        thckGrid(2) = thckGrid(2) + thckGrid(1);
        thckGrid(1) = [];
    end
    
    % (2) %%%%%%%%%%%%% convert Litho model to the correct format %%%%%%%%%%%%%
    maxmodel = 600;    % number of grids in Z
    gridsize = 500;    % thickness of each grid in meter
        
    thckGrid = thckGrid(thckGrid ~= 0);
    rhoGrid = rhoGrid(rhoGrid ~= 0);
    vPGrid = vPGrid(vPGrid ~= 0);
    vSGrid = vSGrid(vSGrid ~= 0);
    
    % compute the depth of each layer
    depthGrid = zeros(length(thckGrid), 1);
    for id = flip(1:length(thckGrid))
        depthGrid(id) = sum(thckGrid(1:id));
    end
    
    depthGrid = round(depthGrid./gridsize);

    if sum(depthGrid) > maxmodel
        % ignore the model deeper than 300m
        depthGrid(end) = maxmodel - sum(depthGrid(1:end-1));
    end
        
    % create the input for Raylee code
    pvs = ones(maxmodel, 1);
    pvp = ones(maxmodel, 1);
    prho = ones(maxmodel, 1);
        
    layerst = 1;
    for i = 1:length(depthGrid)
        pvs(layerst: depthGrid(i)) = pvs(layerst: depthGrid(i)) .* vSGrid(i);
        pvp(layerst: depthGrid(i)) = pvp(layerst: depthGrid(i)) .* vPGrid(i);
        prho(layerst: depthGrid(i)) = prho(layerst: depthGrid(i)) .* rhoGrid(i);
        
        layerst = 1 + depthGrid(i);
    end
        
    % fill to the max. depth
    pvs(layerst: end) = pvs(layerst: end) .* vSGrid(end);
    pvp(layerst: end) = pvp(layerst: end) .* vPGrid(end);
    prho(layerst: end) = prho(layerst: end) .* rhoGrid(end);
    
    % get ADAMA measurement at these points, convert the unit to m/s
    ipavg = pavgs(im,:) .* 1000;
    ipsd = psds(im,:) .* 1000;
    
    try
        % (3) %%%%%%%%%%%%%%%%%%%%%%%%%%% invert ADAMA %%%%%%%%%%%%%%%%%%%%%%%%%%%
        [mA.vsv_update, mA.snsmf_vstot, mA.U, mA.rmserror, mA.chisqurd] = ...
            raylee_invert(pvs.', pvp.', prho.', ipavg, ipsd, maxmodel, gridsize);

        ADAMA1D{im} = mA;  % save the result

    catch ME
        disp(['model ', num2str(im), ' not successful'])
        disp(ME.message)
        notinverse = [notinverse im];
    end
end

% save('./CADAMA_fixed_new.mat', 'ADAMA1D')