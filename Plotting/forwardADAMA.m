%%% ADAMA model is LItho updated using Rayleigh Phase
%%% Now, we use the ADAMA model to forward LP, RG, and LG

addpath(genpath('/Users/sxue3/Documents/Hawaii/Raylee_n_Lovee/raylee_codes-master'))


%% Load models
% use rho from Litho
load('/Users/sxue3/Documents/BayMap_Figures/Data/Litho1D.mat')
% use vs from ADAMA
load('/Users/sxue3/Documents/BayMap_Figures/Data/ADAMAinvert/CADAMA_fixed.mat')

%% forward the model

for im = 1:length(dslon)
    disp(['model: ', num2str(im)])
    
    % (1) %%%%%%%%%%%%% read in ADAMA and Litho models %%%%%%%%%%%%%
    rhoGrid = Litho1D{im}.rho;
    rhoGrid = rhoGrid(1:600);
    
    if isempty(ADAMA{im})   % if ADAMA has no updates
        vSGrid = Litho1D{im}.vs;
        vSGrid = vSGrid(1:600);
        vPGrid = Litho1D{im}.vp;    % NOT SURE: if need to use fixed vp-vs
        vPGrid = vPGrid(1:600);
    else
        ADAMAvs = ADAMA1D{ip}.vsv_update;
        fupdate = find(ADAMAvs(:,1) ~= 0);
        vSGrid = ADAMAvs(fupdate(end),:);
        vPGrid = vSGrid * 1.7;      % fixed vp-vs ratio
        AMupdate = 1;
    end
    
       
    % (2) %%%%%%%%%%%%% create inputs for forward process %%%%%%%%%%%%%
    % construct a grid in solid
    Nn = 600;                   % number of elements in solid (number of grids in Z)
    gridsize = 500;    % thickness of each grid in meter
    h = gridsize*ones(1,Nn);    % grid spacing of mesh (meters)
    % construct a grid in fluid
    Nnf = 0;                  % number of elements in fluid
    hfv = zeros(1,Nnf);          % grid spacing of mesh (meters)
    vpfv = zeros(1,Nnf);
    rhofv = zeros(1,Nnf);

    % forward at the following frequencies:
    period = [5 6 8 10 12 15 20 25 30 35 40];
    fks = 1./period;   
    Nf = length(period);
    % fundamental mode only
    modnv = ones(1,Nf);

    % what type of velocity - phase (0) or group (1)?
    vtypv = zeros(1,Nf);


    % Love Phase Forward


    % Love Group Forward


    % Rayleigh Group Forward


    % save all velocity values to a cell
end