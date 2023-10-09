%%%% Table 2: generate mean velocity and cell numbers for RP, LP, RG, LG
%%%%% Siyu Xue -- May 17. 2023
addpath('/Users/sxue3/Documents/BayMap_Figures/Utils')

% Read in the data
% load Madagascar country boundaries
S1 = load('/Users/sxue3/Documents/BayMap_Figures/Data/AfrCountry.mat').S1;
madgline = polyshape((S1(448).X).',(S1(448).Y).');  % load Madgascar's Big Island

% load the African coast
afcoast = load('/Users/sxue3/Documents/BayMap_Figures/Data/AfricaCoast.mat');
aflon = wrapTo180(afcoast.XY(:, 1));
aflat = afcoast.XY(:,2);
coastline = polyshape(aflon, aflat);

WorkDir = '/Users/sxue3/Documents/BayMap_Figures/Data/LovePData/';
% periods = [5,6,8,10,12,15,20,25,30,35,40];
periods = [5, 25];

for ip = 1:11
    mapname = ['L', num2str(periods(ip))];
    % load the velocity map
    part_xctrs = load([WorkDir, mapname, '/xctrs.mat']).xctrs;
    part_yctrs = load([WorkDir, mapname, '/yctrs.mat']).yctrs;
    part_vels = load([WorkDir, mapname, '/velvals.mat']).velvals;
    
    % Find the in-Africa cells
    
    nummaps = size(part_vels, 1);   % total # of ites
    cellcount = zeros(1, length(1:100:nummaps));    % store the # cells in Africa in each ite
    velavg = cellcount;     % store the avg velocity in each ite
    
    ipoint = 1;
    
    for i = 1:100:nummaps   % loop through every 100 iteration
    
        % x coordinates
        x_thisiter = part_xctrs(i,:);
        x_thisiter = x_thisiter.';
        % y coordinates
        y_thisiter = part_yctrs(i,:);
        y_thisiter = y_thisiter.';
        % cell velocity
        v_thisiter = part_vels(i,:);
        v_thisiter = v_thisiter.';
    
        % check if the cell center is in Africa
        inAfrica = isinterior(coastline, [x_thisiter, y_thisiter]) | isinterior(madgline, [x_thisiter, y_thisiter]) == 1;
    
        % count the # of cells in Africa
        cellcount(ipoint) = sum(inAfrica);
    
        velavg(ipoint) = mean(v_thisiter(inAfrica));
    
        ipoint = ipoint + 1;
    end
    
    disp('-------------------------------------')
    disp(mapname)
    disp(['Avg Vel is: ', num2str(mean(velavg))])
    disp(['Avg Cell Count is: ', num2str(mean(cellcount))])
    disp('-------------------------------------')

end