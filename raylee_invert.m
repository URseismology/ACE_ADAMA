%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% PROGRAM:
% raylee_invert.m
%
% PROGRAMMERS:
% Matt Haney and Victor Tsai
%
% Last revision date:
% 26 April 2017
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This code is distributed as part of the source-code package 
%                   raylee_inversion_codes 
% that accompanies Haney and Tsai (2017). The package can be downloaded 
% from the Geophysics source-code archive at 
%                   http://software.seg.org/2017/0003/index.html
% Use of this code is subject to acceptance of the terms and conditions
% that can be found at http://software.seg.org/disclaimer.txt 
% Copyright 2017 by The Society of Exploration Geophysicists (SEG)
% Reference:
% Haney, M. M., Tsai, V. C. (2017) Perturbational and nonperturbational 
% inversion of Rayleigh-wave velocities, Geophysics, 82(3), F15-F28.
% doi: 10.1190/geo2016-0397.1
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Program raylee_invert is a Matlab script to invert any collection of 
% fundamental mode/higher mode Rayleigh wave group or phase 
% velocities measured at a set of frequencies for a shear wave velocity 
% depth model.
%
% The program can invert for shear velocity assuming the Vp/Vs ratio is 
% fixed in the subsurface or assuming that the original Vp is unchanged
% during the inversion.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Input
%
% 1. Input parameters are read from the file "input_params.txt"
%
% 2. Measured (group or phase) velocities are read from 
%    the file "velocity_values.txt"
%
% 3. The frequencies at which the velocities are measured are read 
%    from the file "frequency_values.txt"
%
% 4. Error bars on the measured (group or phase) velocities are read from 
%    the file "velocity_values_errs.txt"
%
% 5. The types of velocities measured (group or phase) are read from 
%    the file "vtype_values.txt"
%
% 6. The mode numbers of the measurements (group or phase) are read from 
%    the file "mode_values.txt"
%
% 7. The finite-element grid used for modeling and inversion is read 
%    from the files "grid_values_solid.txt" and "grid_values_fluid.txt"
%
% 8. The initial model of the subsurface material properties 
%    (Vp, Vs, and density) for the inversion is read from the files
%    "vp_init.txt", "vs_init.txt", and "rho_init.txt". 
%
% 9. The subsurface material properties (Vp and density) for a water 
%    layer above the solid are read from the files "vpf.txt" and 
%    "rhof.txt". 
%    
% Output
%
% 1. The Vs model updates are available in the matrix vsv_update. This 
%    matrix is size (nupdat x Nn) where nupdat is the number of iterations 
%    executed before the stopping criterion is met and Nn is the number 
%    of nodes for the finite-element grid. The first update is in row 1 
%    and the final update is in row nupdat. 
%
% 2. The sensitivity kernel matrix for the last iteration is the 
%    matrix snsmf_vstot. This matrix has size (Nn x Nf), where Nf is the 
%    number of frequencies at which measurements exist.
%
% 3. The computed group or phase velocity for the final update is the 
%    vector U, size (1 x Nf).
%
% 4. The RMS error for the initial guess and all the updates is stored 
%    in the vector rmserror, size (1 x (nupdat+1)).
%
% 5. The Chi-squared error for the initial guess and all updates is stored 
%    in the vector chisqurd, size (1 x (nupdat+1)).
%

function [vsv_update, snsmf_vstot, U, rmserror, chisqurd] = ...
    raylee_invert(vsv, vpv, rhov, U_data, U_data_errs, Nn, gridsize)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this is a script
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% begin input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% vsv = vector of initial Vs model (Litho)
% vpv = vector of initial Vp model (Litho)
% rhov = vector of initial density model (Litho)
% U_data = ADAMA measurement (period short to long)
% U_data_errs = error bars of ADAMA velocity measurement
% Nn = number of elements in solid

% inversion parameters
lsmth = gridsize; % Siyu: not sure      % model smoothness scale (m)
msigmaf = 2;      % Siyu: not sure      % model standard deviation factor 
pratioflag = 1;                         % 1: fix vp - vs ratio

if pratioflag
    vpv = 0.9409e3 + vsv.*(2.0947 + vsv.*(-0.8206e-3 + vsv.*(0.2683e-6 - 0.0251e-9*vsv)));
    rhov = vpv.*(1.6612 + vpv.*(-0.4721e-3 + vpv.*(0.0671e-6 + vpv.*(-0.0043e-9 + 0.000106e-12*vpv))));
end
                                     
% stopping criteria
nupds = 16;   % Siyu: not sure     % max number of updates (iterations)

% chi squared window
chilo = 0.3;   % Siyu: not sure
chihi = 0.5;      % Siyu: not sure

% load grid
h = gridsize*ones(1,Nn);             % grid spacing of mesh (meters)

% load grid in fluid
Nnf = 0;       % number of elements for fluid
hfv = 100*ones(1,Nnf);

% load frequencies
pds = [5 6 8 10 12 15 20 25 30 35 40];
fks = 1./pds;   % get the frequencies

Nf = 11;                        % number of measurements
modn = ones(1, Nf);             % data vector of mode numbers
vflg = zeros(1, Nf);            % data vector of velocity types (0=phase)

% load Vp model in fluid
vpfv = zeros(1,Nnf);

% load density model in fluid
rhofv = zeros(1,Nnf);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% end input 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% initialize velocity update vector
vsv_update = zeros(nupds,Nn);

% make a vector of depths at nodes by a running sum of the grid spacings
hs(1) = 0;
for ii=2:length(h)
    hs(ii) = sum(h(1:(ii-1)));
end

% make a vector of depths at center of elements by a running sum 
hss(1) = h(1)/2;
for ii=2:length(h)
    hss(ii) = sum(h(1:(ii-1))) + h(ii)/2;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% check whether input parameters are physically possible
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% density greater than zero
if (sum(rhov <= 0) > 0)
    error('Negative density values exist in initial guess');
else
end
% shear velocity greater than zero
if (sum(vsv <= 0) > 0)
    error('Negative shear velocity values exist in initial guess');
else
end
% poisson's ratio between two bounds
pratio = (vpv.^2 - 2*(vsv.^2))./(2*(vpv.^2 - vsv.^2));
if ((sum(pratio <= -1) > 0) || (sum(pratio >= 0.5) > 0))
    error('Impossible Poisson ratio values exist in initial guess');
else
end
% density greater than zero in fluid
if (sum(rhofv <= 0) > 0)
    error('Negative density values exist in initial guess');
else
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% prepare for initial inversion step
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% compute sensitivity kernel using initial guess
[U, snsmf_vstot] = raylee_sensitivity(Nn,vsv,vpv,...
                    rhov,fks,h,modn,vflg,Nnf,vpfv,rhofv,hfv,pratioflag);

% find the measurements for which both data and model are not NaN
[Ur, U_datar, fksr, fksri, modnr, vflgr, snsmf_vstotr] = ...
    check_nans(U, U_data, fks, modn, vflg, snsmf_vstot);
Nfr = length(fksr);

% save the S-wave velocity guess and the resulting data
vsv_guess = vsv;
U_guess = Ur;
fksr_guess = fksr;

% calculate the a priori model covariance matrix and inverse square root
msigma = mean(U_data_errs(fksri))*msigmaf;
mcm = (msigma^2)*exp(-abs(repmat(hs,Nn,1)-repmat(hs',1,Nn))/lsmth);
mcmisr = sqrtm(inv(mcm));

% calculate the a priori data covariance matrix and inverse square root
dcm = diag(U_data_errs(fksri).^2);
dcmisr = diag(1./U_data_errs(fksri));

% rms error of the initial guess
rmserror(1) = sqrt(mean(((U_guess-U_datar)./1).^2));
chisqurd(1) = (U_guess-U_datar)*dcmisr*dcmisr*transpose(U_guess-U_datar);
Nfrv(1) = Nfr;

% check to see if initial guess has chi^2 less than 1
if ((chisqurd(1)/Nfr) < chilo)
    error('Initial model fits data to less than 1 chi-squared');
elseif ((chisqurd(1)/Nfr) < chihi)
    error('Initial model fits data within acceptable chi-squared window');
else    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% invert using damped least squares method of Tarantola and Valette (1982)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% initial damped linear inversion
dvs = linvers(U_datar,Ur,snsmf_vstotr,mcmisr,dcmisr,Nn,vsv,vsv_guess);

% add to the initial model
vsv = dvs' + vsv_guess;
if (pratioflag == 1)
    vpv = 0.9409e3 + vsv.*(2.0947 + vsv.*(-0.8206e-3 + vsv.*(0.2683e-6 - 0.0251e-9*vsv)));
    rhov = vpv.*(1.6612 + vpv.*(-0.4721e-3 + vpv.*(0.0671e-6 + vpv.*(-0.0043e-9 + 0.000106e-12*vpv))));
end

% compute new sensitivity kernel
[U, snsmf_vstot] = raylee_sensitivity(Nn,vsv,vpv,...
                    rhov,fksr,h,modnr,vflgr,...
                    Nnf,vpfv,rhofv,hfv,pratioflag);

% find NaNs
[U, U_datar, fksr, fksri, modnr, vflgr, snsmf_vstot] = ...
    check_nans(U, U_datar, fksr, modnr, vflgr, snsmf_vstot);

% if number of NaNs changed, recompute data and model covariances
if (length(fksr) ~= Nfr)
    Nfr = length(fksr);
    msigma = mean(U_data_errs(fksri))*msigmaf;
    mcm = (msigma^2)*exp(-abs(repmat(hs,Nn,1)-repmat(hs',1,Nn))/lsmth);
    mcmisr = sqrtm(inv(mcm));
    dcm = diag(U_data_errs(fksri).^2);
    dcmisr = diag(1./U_data_errs(fksri));
else
end
                
% compute RMS error and chi-squared
rmserrorp = sqrt(mean(((U-U_datar)./1).^2));
chisqurdp = (U-U_datar)*dcmisr*dcmisr*transpose(U-U_datar); 
    
% a reduced line search if chi^2 of update is not lower
nreds = 0;
while ((chisqurdp >= chisqurd(1) && nreds < nupds) || ...
        ((chisqurdp/Nfr) < 1 && nreds < nupds))
        
    nreds = nreds + 1;
        
    % reduce step by a factor of 2, and add it in
    dvs = dvs/2;
    vsv = vsv_guess + dvs';
    if (pratioflag == 1)
        vpv = 0.9409e3 + vsv.*(2.0947 + vsv.*(-0.8206e-3 + vsv.*(0.2683e-6 - 0.0251e-9*vsv)));
        rhov = vpv.*(1.6612 + vpv.*(-0.4721e-3 + vpv.*(0.0671e-6 + vpv.*(-0.0043e-9 + 0.000106e-12*vpv))));
    end
    
    % call the sensitivity function to compute U
    [U, snsmf_vstot] = raylee_sensitivity(Nn,vsv,vpv,...
                        rhov,fksr,h,modnr,vflgr,...
                        Nnf,vpfv,rhofv,hfv,pratioflag);
                    
    % check for NaNs
    [U, U_datar, fksr, fksri, modnr, vflgr, snsmf_vstot] = ...
        check_nans(U, U_datar, fksr, modnr, vflgr, snsmf_vstot);
    
    % if number of NaNs changed recompute data and model covariances
    if (length(fksr) ~= Nfr)
        Nfr = length(fksr);
        msigma = mean(U_data_errs(fksri))*msigmaf;
        mcm = (msigma^2)*exp(-abs(repmat(hs,Nn,1)-repmat(hs',1,Nn))/lsmth);
        mcmisr = sqrtm(inv(mcm));
        dcm = diag(U_data_errs(fksri).^2);
        dcmisr = diag(1./U_data_errs(fksri));
    else
    end

    % the rms of this potential update
    rmserrorp = sqrt(mean(((U-U_datar)./1).^2));
    % the chi^2 of this potential update
    chisqurdp = (U-U_datar)*dcmisr*dcmisr*transpose(U-U_datar);
    
end
    
% shear velocity must be greater than zero
if (sum(vsv <= 0) > 0)
    error('Negative shear velocity values encountered in inversion');
else
end
% poisson's ratio between two bounds
pratio = (vpv.^2 - 2*(vsv.^2))./(2*(vpv.^2 - vsv.^2));
if ((sum(pratio <= -1) > 0) || (sum(pratio >= 0.5) > 0))
    error('Impossible Poisson ratio values encountered in inversion');
else
end
    
% the updated model, print number of update to screen
nupdat = 1;
vsv_update(nupdat,:) = vsv; 
    
% the rms of this update
rmserror(nupdat+1) = rmserrorp;
% the chi^2 of this update
chisqurd(nupdat+1) = chisqurdp;

end

% now full modeling
[U, snsmf_vstot] = raylee_sensitivity(Nn,vsv,vpv,...
                    rhov,fks,h,modn,vflg,...
                    Nnf,vpfv,rhofv,hfv,pratioflag);
                
% check for NaNs  
[Ur, U_datar, fksr, fksri, modnr, vflgr, snsmf_vstotr] = ...
    check_nans(U, U_data, fks, modn, vflg, snsmf_vstot);

% if number of NaNs changed recompute data and model covariances
if (length(fksr) ~= Nfr)
    Nfr = length(fksr);
    msigma = mean(U_data_errs(fksri))*msigmaf;
    mcm = (msigma^2)*exp(-abs(repmat(hs,Nn,1)-repmat(hs',1,Nn))/lsmth);
    mcmisr = sqrtm(inv(mcm));
    dcm = diag(U_data_errs(fksri).^2);
    dcmisr = diag(1./U_data_errs(fksri));    
else
end

% compute RMS and chi-squared
rmserrorp = sqrt(mean(((Ur-U_datar)./1).^2));
chisqurdp = (Ur-U_datar)*dcmisr*dcmisr*transpose(Ur-U_datar);

% the rms of this update
rmserror(nupdat+1) = rmserrorp;
% the chi^2 of this update
chisqurd(nupdat+1) = chisqurdp;
Nfrv(nupdat+1) = Nfr;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% now an iterative loop, updating the initial guess
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% while the stopping criterion and the maximum 
% allowed number of iterations has not been met, continue updating
while ((chisqurdp/Nfr) > chihi && nupdat < nupds ) 

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % invert again as in Tarantola and Valette (1982)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
    
    % linear inverse
    dvs = linvers(U_datar,Ur,snsmf_vstotr,mcmisr,dcmisr,Nn,vsv,vsv_guess);
    
    % add to the initial model
    vsv = dvs' + vsv_guess;
    % if fixed vpvs ratio, adjust vp model
    if (pratioflag == 1)
        vpv = 0.9409e3 + vsv.*(2.0947 + vsv.*(-0.8206e-3 + vsv.*(0.2683e-6 - 0.0251e-9*vsv)));
        rhov = vpv.*(1.6612 + vpv.*(-0.4721e-3 + vpv.*(0.0671e-6 + vpv.*(-0.0043e-9 + 0.000106e-12*vpv))));
    else
    end
    
    % call the sensitivity function to model
    [U, snsmf_vstot] = raylee_sensitivity(Nn,vsv,vpv,...
                    rhov,fksr,h,modnr,vflgr,...
                    Nnf,vpfv,rhofv,hfv,pratioflag);

    % check for NaNs
    [U, U_datar, fksr, fksri, modnr, vflgr, snsmf_vstot] = ...
        check_nans(U, U_datar, fksr, modnr, vflgr, snsmf_vstot);
    
    % if number of data changed, recompute data and model covariances
    if (length(fksr) ~= Nfr)
        Nfr = length(fksr);
        msigma = mean(U_data_errs(fksri))*msigmaf;
        mcm = (msigma^2)*exp(-abs(repmat(hs,Nn,1)-repmat(hs',1,Nn))/lsmth);
        mcmisr = sqrtm(inv(mcm));
        dcm = diag(U_data_errs(fksri).^2);
        dcmisr = diag(1./U_data_errs(fksri));        
    else
    end                
                    
    % compute rms and chi of this potential update
    rmserrorp = sqrt(mean(((U-U_datar)./1).^2));
    chisqurdp = (U-U_datar)*dcmisr*dcmisr*transpose(U-U_datar); 
    
    % a reduced line search if chi^2 of update is not lower
    nreds = 0;
    % the gradient - difference between the current update and previous
    dvs = (vsv' - transpose(vsv_update(nupdat,:)));
    
    while ((chisqurdp >= 1.01*chisqurd(nupdat+1) && nreds < nupds) || ...
            ((chisqurdp/Nfr) < chilo && nreds < nupds))
        
        nreds = nreds + 1;
        
        % reduce step by a factor of 2, and add it in
        dvs = dvs/2;
        vsv = vsv_update(nupdat,:) + dvs';
        % if vpvs ratio fixed, adjust vp model
        if (pratioflag == 1)
            vpv = 0.9409e3 + vsv.*(2.0947 + vsv.*(-0.8206e-3 + vsv.*(0.2683e-6 - 0.0251e-9*vsv)));
            rhov = vpv.*(1.6612 + vpv.*(-0.4721e-3 + vpv.*(0.0671e-6 + vpv.*(-0.0043e-9 + 0.000106e-12*vpv))));
        end
    
        % call the sensitivity function to compute U
        [U, snsmf_vstot] = raylee_sensitivity(Nn,vsv,vpv,...
                            rhov,fksr,h,modnr,vflgr,...
                            Nnf,vpfv,rhofv,hfv,pratioflag);
        
        % check for NaNs
        [U, U_datar, fksr, fksri, modnr, vflgr, snsmf_vstot] = ...
            check_nans(U, U_datar, fksr, modnr, vflgr, snsmf_vstot);
        %Nfr = length(fksr);
        
        % if number of data changed, adjust model and data covariances
        if (length(fksr) ~= Nfr)
            Nfr = length(fksr);
            msigma = mean(U_data_errs(fksri))*msigmaf;
            mcm = (msigma^2)*exp(-abs(repmat(hs,Nn,1)-repmat(hs',1,Nn))/lsmth);
            mcmisr = sqrtm(inv(mcm));
            dcm = diag(U_data_errs(fksri).^2);
            dcmisr = diag(1./U_data_errs(fksri));           
        else
        end                    
        
    % the rms of this potential update
    rmserrorp = sqrt(mean(((U-U_datar)./1).^2));
    % the chi^2 of this potential update
    chisqurdp = (U-U_datar)*dcmisr*dcmisr*transpose(U-U_datar);
        
    end
    
    % shear velocity must be greater than zero
    if (sum(vsv <= 0) > 0)
        error('Negative shear velocity values encountered in inversion');
    else
    end
    % poisson's ratio between two bounds
    pratio = (vpv.^2 - 2*(vsv.^2))./(2*(vpv.^2 - vsv.^2));
    if ((sum(pratio <= -1) > 0) || (sum(pratio >= 0.5) > 0))
        error('Impossible Poisson ratio values encountered in inversion');
    else
    end
    
    % the next updated model, print number of update to screen
    nupdat = nupdat + 1;
    vsv_update(nupdat,:) = vsv; 
    
    % the rms of this update
    rmserror(nupdat+1) = rmserrorp;
    % the chi^2 of this update
    chisqurd(nupdat+1) = chisqurdp;
    
    % now full modeling
    [U, snsmf_vstot] = raylee_sensitivity(Nn,vsv,vpv,...
                    rhov,fks,h,modn,vflg,...
                    Nnf,vpfv,rhofv,hfv,pratioflag);
                
    % check for NaNs   
    [Ur, U_datar, fksr, fksri, modnr, vflgr, snsmf_vstotr] = ...
        check_nans(U, U_data, fks, modn, vflg, snsmf_vstot);
    %Nfr = length(fksr);

    % if number of data changed, recompute data and model covariances
    if (length(fksr) ~= Nfr)
        Nfr = length(fksr);
        msigma = mean(U_data_errs(fksri))*msigmaf;
        mcm = (msigma^2)*exp(-abs(repmat(hs,Nn,1)-repmat(hs',1,Nn))/lsmth);
        mcmisr = sqrtm(inv(mcm));
        dcm = diag(U_data_errs(fksri).^2);
        dcmisr = diag(1./U_data_errs(fksri));    
    else
    end

    % compute rms and chi^2
    rmserrorp = sqrt(mean(((Ur-U_datar)./1).^2));
    chisqurdp = (Ur-U_datar)*dcmisr*dcmisr*transpose(Ur-U_datar);

    % the rms of this update
    rmserror(nupdat+1) = rmserrorp;
    % the chi^2 of this update
    chisqurd(nupdat+1) = chisqurdp;
    Nfrv(nupdat+1) = Nfr;
                
end

% end the timer

sprintf('%d of %d measurements used',Nfr,Nf-sum(isnan(U_data)))

if ((chisqurd(nupdat+1)/Nfr) > chihi)
    sprintf('WARNING: Inversion did not converge to stopping criterion and underfitted data. Increase number of updates.')
else
end

if ((chisqurd(nupdat+1)/Nfr) < chilo)
    sprintf('WARNING: Inversion did not converge to stopping criterion and overfitted data. Increase number of reduction steps.')
else
end






