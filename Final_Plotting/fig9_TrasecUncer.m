%%%%%%%% Compute the uncertanity for ADAMA velocity-depth updates
%%% generate data for Fig. 9
%%% Siy Xue - Apr 6 2023
addpath(genpath('/Users/sxue3/Documents/Hawaii/Raylee_n_Lovee/raylee_codes-master'))
addpath('/Users/sxue3/Documents/BayMap_Figures/Utils')

%% Simulate and compute inversion uncertainty

% index of downsampled points in the transec
tsi = [1068 1365 2297 2317 530 1034 291 1255 661 1703 1925 1184 370 1707 145 ...
    2071 731 1031 1861 1642 602 1195 232 973 382 1720 1799 1325 1027 1857 ...
    79 2201 415 1752 2275 2340 1392 795 1539];

% Load the ADAMA crust updates
load('/Users/sxue3/Documents/BayMap_Figures/Data/ADAMAinvert/ADAMA_1D_fixedpoisson.mat')

% Load the ADAMA downsample locations
dsdots = load('/Users/sxue3/Documents/BayMap_Figures/Data/AfrDownSampleDots.mat');
dslat = dsdots.dslat;
dslon = dsdots.dslon;

% Load the Litho reference models
load('/Users/sxue3/Documents/BayMap_Figures/Data/Litho1D.mat')

% Load the velocity maps sample locations
load('/Users/sxue3/Documents/BayMap_Figures/Data/RayleighPData/plotpara.mat');
longrid = longrid(1,:);
latgrid = latgrid(:,1);

% Get the sample points in map
afrmask = zeros(2380, 2);
for i = 1:length(dslat)
    xi = find(latgrid <= dslat(i), 1, 'last');
    yi = find(longrid <= dslon(i), 1, 'last');

    afrmask(i, :) = [xi, yi];
end

% loop through all points
for id = 1:2381
    if ~isempty(ADAMA1D{id})
%         % simulate and save non-ADAMA updates
%         disp(['Sampling dot: ' num2str(id)])
%         getUncer(id, Litho1D, afrmask);
        continue
    else
        % simulate and save ADAMA updated models
        disp(['Sampling dot: ' num2str(id)])
        getUncer(id, Litho1D, afrmask);
    end
end


%% plot the simulation result of a sampled location (if needed)
hold on

ADAMAidot = [];
for im = 1:100

    if ~isempty(ADAMAid{im})
        ADAMAvs = ADAMAid{im}.vsv_update;
        fupdate = find(ADAMAvs(:,1) ~= 0);
        plot(ADAMAvs(fupdate(end),:)./1000, 'r')
        ADAMAidot = [ADAMAidot; ADAMAvs(fupdate(end),1:120)./1000];
    end 
end

% Compute the SD at each 500m depth
velsd = std(ADAMAidot, 0, 1);
velmean = mean(ADAMAidot, 1);

plot(velmean, 'k', 'LineWidth',2)
plot(velmean+velsd, 'b', 'LineWidth',1.5)
plot(velmean-velsd, 'b', 'LineWidth',1.5)
xlim([1,120])
xticks([0 20 40 60 80 100])
xticklabels({'0','10','20','30','40','50'})  % grid size is 500m
xlabel('Depth [km]')
ylabel('Vel [km/s]')

%% assemble all uncertanities together (if needed)



%% Function for simulation
function getUncer(idot, Litho1D, afrmask)

    type = 'P';
    Periods = {'R5', 'R6', 'R8', 'R10', 'R12', 'R15', 'R20', 'R25', 'R30', 'R35', 'R40'};
    
    for ip = 1:length(Periods)
        period = Periods{ip};
        
        pavgpath = ['/Users/sxue3/Documents/BayMap_Figures/Data/RayAvgMap/', period, '_', type, '_maps.mat'];
        pavg = load(pavgpath).avgmap;
        psd = load(pavgpath).sdevmap;
        
        ipavg = pavg(afrmask(idot, 1), afrmask(idot, 2));
        ipsd = psd(afrmask(idot, 1), afrmask(idot, 2));
    
        pavgs(:,ip) = ipavg;
        psds(:,ip) = ipsd;
    end
    
    
    %%% Shift the mean curve up and down randomly
    
    velrefs = pavgs;
    temp = linspace(pavgs(1), pavgs(4),4);
    % velrefs(1:4) = temp;
    
    x = 1:11;
    y = velrefs;
    xx = 1:.1:11;
    yyvel = interp1(x,y,xx, 'spline');
    plot(x,y,'o',xx,yyvel)
    % 
    % yysd = interp1(x,psds,xx, 'spline');
    % plot(x,psds,'o',xx,yysd)
    
    % Sampling the dispersion curves
    dcsamples = zeros(11, 100);
    for ic = 1:100 
        temp = randi([1 11],1);
        dcsamples(:, ic) = randn() * psds(temp) + pavgs;
    end
    
    % Plot to check the sampling
    % hold on
    % for iii = 1:100
    %     plot(x, dcsamples(:, iii), 'r')
    %     pause(0.1)
    % end
    % 
    % errorbar(pavgs, psds, 'k.');    % ADAMA measurements w. uncert.
    
    %%% Raylee Inversion
    
    ADAMAid = cell(100, 1);
    notinverse = [];
     
    maxmodel = 600;    % number of grids in Z
    gridsize = 500;    % thickness of each grid in meter
    
    pvs = Litho1D{idot}.pvs;
    pvp = Litho1D{idot}.pvp;
    prho = Litho1D{idot}.prho;
    
    pvs = pvs(1:maxmodel);
    pvp = pvp(1:maxmodel);
    prho = prho(1:maxmodel);
    
    for im = 1:100
        disp(['model: ', num2str(im)])
        
        % get ADAMA measurement at these points, convert the unit to m/s
        ipavg = dcsamples(:,im)' .* 1000;
        ipsd = psds * 1000;
        
        try
            % (3) %%%%%%%%%%%%%%%%%%%%%%%%%%% invert ADAMA %%%%%%%%%%%%%%%%%%%%%%%%%%%
            [mA.vsv_update, mA.snsmf_vstot, mA.U, mA.rmserror, mA.chisqurd] = ...
                raylee_invert(pvs.', pvp.', prho.', ipavg, ipsd, maxmodel, gridsize);
    
            ADAMAid{im} = mA;  % save the result
    
        catch ME
            disp(['model ', num2str(im), ' not successful'])
            disp(ME.message)
            notinverse = [notinverse im];
        end
    end
    
    save(['/Users/sxue3/Documents/BayMap_Figures/Data/TransecUncer/', num2str(idot), '.mat'], 'ADAMAid', 'dcsamples', 'pavgs', 'psds')

end


