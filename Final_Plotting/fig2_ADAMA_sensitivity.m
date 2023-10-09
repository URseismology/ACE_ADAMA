%%%% Fig 2: plot the depth sensitivity for Love and Rayleigh

% For Yuri make rayleigh sensitivity...
% Dr. O

% this is a script'

% model from XD.RUNG in tanzania 5 
%  0.019000 0   1.070000   0.010700
% 12.185000 0   3.573320   0.035733
% 14.491000 0   3.711490   0.037115
% 14.583000 0   3.930270   0.039303
%  0.000000 0   4.596920   0.045969

%%%% modified: Siyu Xue -- Feb 10. 2023
%%% Note: this script plots for both Fig2 & Fig 4

%% --- specify target model here -- some standard LITHO model
clear all
clc
addpath(genpath('/Users/sxue3/Documents/Hawaii/Raylee_n_Lovee'))
addpath('/Users/sxue3/Documents/BayMap_Figures/Utils')

% construct a grid in solid
Nn = 240;                   % number of elements in solid
Nn = 1000;
h = 250*ones(1,Nn);             % grid spacing of mesh (meters)
% construct a grid in fluid
Nnf = 0;                  % number of elements in fluid
hfv = 100*ones(1,Nnf);          % grid spacing of mesh (meters)

% construct vector of frequencies
Nf = 56;                   % number of measurements
fmin = 0.10;                    % minimum frequency (Hz)
df = 0.01;                      % frequency spacing (Hz)
fks = fmin+(df*[0:(Nf-1)]);     % vector of frequencies (Hz)
                                % these are the frequencies at which 
                                % the velocities are measured
                                
% ADAMA1 measurement frequency (see Figure 8)
Nf = 56;
fmin = 25e-3; % 25 mHz
fmax = 200e-3; % 200 mHz
fks = linspace(fmin, fmax, Nf);
                            
% -- model 2 modes, 1 fundamental and 1st higher mode 
%Nf = 2*Nf;
% vector of mode numbers
%modnv = [ones(1,Nf/2) 2*ones(1,Nf/2)];
% try to model over same frequencies
%fks = [fks fks]; 
% what type of velocity - phase (0) or group (1)?
%vtypv = zeros(1,Nf);

% using only fundamental modes in ADAMA - ambient noise
% specify vector of mode number for each frequency measurements
modn = ones(1,Nf);  % mode number here -- fundamental
vflg = zeros(1,Nf); % code for phase or group velocity --use 0 for phase vel
pratioflag = 1; % fixed vp vs ratio


% --
% medium parameters
vpvsr = 1.7321;                 % Vp/Vs ratio
gardc = 309.6;                 % constant in Gardner relation
powr = 0.25;                  % exponent in Gardner relation

% make a three layered model
layrth1 = 12; % thickness in elements, layrth1*h = thickness in meters
layrth2 = 15; % thickness in elements, layrth2*h = thickness in meters
layrth3 = 15;
% the true model in the solid
vslay1 = 3573.320 ; vplay1 = vslay1*vpvsr; rholay1 = gardc*(vplay1^powr);
vslay2 = 3711.490; vplay2 = vslay2*vpvsr; rholay2 = gardc*(vplay2^powr);
vslay3 = 3930.270; vplay3 = vslay3*vpvsr; rholay3 = gardc*(vplay3^powr);
vslay4 = 4596.920; vplay4 = vslay4*vpvsr; rholay4 = gardc*(vplay4^powr);
vpv = [vplay1*ones(1,layrth1) vplay2*ones(1,layrth2) ...
    vplay3*ones(1,layrth3) ...
       vplay4*ones(1,(Nn-(layrth1+layrth2+layrth3)))];
vsv = [vslay1*ones(1,layrth1) vslay2*ones(1,layrth2) ...
    vslay3*ones(1,layrth3) ...
       vslay4*ones(1,(Nn-(layrth1+layrth2+layrth3)))];
rhov = [rholay1*ones(1,layrth1) rholay2*ones(1,layrth2) ...
    rholay3*ones(1,layrth3) ...
        rholay4*ones(1,(Nn-(layrth1+layrth2+layrth3)))];

% make the fluid part of the model    
vplay4 = 1500; rholay4 = 1000;
% the true model in the fluid
vpfv = vplay4*ones(1,Nnf);
rhofv = rholay4*ones(1,Nnf);    

% %%%%%%%% -- done specifying model

% compute sensitivity kernel using initial guess
[U, ral_snsmf_vstot] = raylee_sensitivity(Nn,vsv,vpv,...
                    rhov,fks,h,modn,vflg,Nnf,vpfv,rhofv,hfv,pratioflag);

[~, love_snsmf_vstot] = lovee_sensitivity(Nn,vsv,rhov,fks,h,modn,vflg);

% make a vector of depths at center of elements by a running sum 
hss(1) = h(1)/2;
for ii=2:length(h)
    hss(ii) = sum(h(1:(ii-1))) + h(ii)/2;
end

% love sensitivity interpolated unto a nice grid
% interpolate onto regular grid
snsmf_vstoti = zeros(length([0:min(h):sum(h)]),Nf);
for ii=1:Nf
snsmf_vstoti(:,ii) = interp1(hss,love_snsmf_vstot(:,ii),[0:min(h):sum(h)],'linear');
end

%% plot sensitivity of reference Litho model
% sensitivity kernel of final update
clf
%figure
fsize = 14;

ax1 = axes();
[qq zz] = size(ral_snsmf_vstot);
smax = round(max(max(ral_snsmf_vstot/min(h)))*10000)/10;
smin = round(min(min(ral_snsmf_vstot/min(h)))*10000)/10;

% % find 56 and 109
% % subplot(1,2,1)
% imagesc(ax1, 1./fks(1:56),0.001*hss,ral_snsmf_vstot(:,1:56)./repmat(.001*h',1,56)); 
% colormap('hot');
% %caxis([0 0.1])
% axis([5 40 0 50])
% hh = colorbar('EastOutside','FontSize',fsize,...
%     'FontWeight','bold','Ytick',[smin:((smax-smin)/5):smax],'Ylim',[smin smax]);
% label = sprintf(' Sensitivity (km^{-1}) ');
% set(get(hh,'YLabel'),'String',label,'FontSize',fsize,'FontWeight','bold')
% set(gca,'Fontsize',fsize,'FontWeight','bold');
% xlabel(' Period (seconds) '); ylabel(' Depth (km) '); 
% title(' V_{S} Kernel Rayleigh Phase Velocity')
% 
% hold on
% xline(28.6, 'w');
% % arrows
% annotation('doublearrow',[0.81 0.125],[0.2 0.2],'Color', 'w', 'LineWidth',2)
% annotation('textbox',[.43 .17 .15 .1], 'String', 'ADAMA','EdgeColor','none', 'FontSize', 14, 'Color', 'w')
% annotation('arrow',[0.59 0.81],[0.3 0.3],'Color', 'w', 'LineWidth',2)
% annotation('textbox',[.63 .27 .2 .1], 'String', 'Litho1.0','EdgeColor','none', 'FontSize', 14, 'Color', 'w')

% -------------
% subplot(1,2,2)
[qq zz] = size(snsmf_vstoti);
smax = round(max(max(snsmf_vstoti/min(h)))*10000)/10;
smin = round(min(min(snsmf_vstoti/min(h)))*10000)/10;
%imagesc(1./fks,hss/1000,1000*snsmf_vstoti(:,1:81)/min(h)); colormap('jet'); caxis([smin smax]); hold on
imagesc(1./fks,hss/1000,1000*snsmf_vstoti/min(h)); colormap('hot'); caxis([smin smax]); hold on
axis([5 40 0 50])
hh = colorbar('EastOutside','FontSize',fsize,...
    'FontWeight','bold','Ytick',[smin:((smax-smin)/5):smax],'Ylim',[smin smax]);
label = sprintf(' Sensitivity (km^{-1}) ');
set(get(hh,'YLabel'),'String',label,'FontSize',fsize,'FontWeight','bold')
set(gca,'Fontsize',fsize,'FontWeight','bold');
xlabel(' Period (seconds) '); ylabel(' Depth (km) '); 
title(' V_{S} Kernel Love Phase Velocity ')

hold on
xline(33.3, 'w');
% arrows
annotation('doublearrow',[0.81 0.125],[0.2 0.2],'Color', 'w', 'LineWidth',2)
annotation('textbox',[.43 .17 .15 .1], 'String', 'ADAMA','EdgeColor','none', 'FontSize', 14, 'Color', 'w')
annotation('arrow',[0.68 0.81],[0.3 0.3],'Color', 'w', 'LineWidth',2)
annotation('textbox',[.685 .27 .2 .1], 'String', 'Litho1.0','EdgeColor','none', 'FontSize', 14, 'Color', 'w')


% set the size of the plot
x0=10;
y0=10;
fwidth=500;
fheight=500;
set(gcf,'position',[x0,y0,fwidth,fheight]);

% save the figure
figpath = '/Users/sxue3/Documents/BayMap_Figures/fig/fig2_LoveSensitivity.pdf';
set(gcf,'Units','Inches');
pos = get(gcf,'Position');
set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
print(gcf,figpath,'-painters', '-dpdf','-r0');