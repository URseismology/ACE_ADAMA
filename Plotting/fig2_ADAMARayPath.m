%%%% Fig 2: ADAMA & Litho raypath maps
%%%% Siyu Xue -- Jan 10. 2023
addpath('/Users/sxue3/Documents/BayMap_Figures/Utils')
addpath('/Users/sxue3/Documents/Hawaii/spheretri-master/');
addpath('/Users/sxue3/Documents/ADAMA_Figures/m_map/');

%% 1. path file 
inPathFile = '/Users/sxue3/Documents/ADAMA_Figures/data/ADAMA_staconns.csv';

% load africa coverage  ...
fname = ['/Users/sxue3/Documents/BayMap_Figures/Data/AfricaCoast.mat'];
afcoast = load(fname);
aflon = wrapTo180(afcoast.XY(:,1));
aflat = afcoast.XY(:,2);
coastline = polyshape(aflon, aflat);

S1 = load('/Users/sxue3/Documents/BayMap_Figures/Data/AfrCountry.mat').S1;
Malon = (S1(448).X).'; Malon = Malon(1:10:end);
Malat = (S1(448).Y).'; Malat = Malat(1:10:end);
madgline = polyshape((S1(448).X).',(S1(448).Y).');  % load Madgascar's Big Island

dat = readtable(inPathFile);
lon1 = dat.lon1;
lat1 = dat.lat1;
lon2 = dat.lon2;
lat2 = dat.lat2;
index = 1:length(lon1);

ploc = [index' lon1 lat1 lon2 lat2];

%2.  discretization of final raycoverage
delP = 0.025;   % path delta
%3. discretization of final map
delDeg = 0.25;  % map delta

%4. re-weigh based on max distance
maxDeg = 10; %maximum distance to rescale -- 30 degrees

%% If need to compile the ADAMA raypath ... 
% %% 5. confirm face triangulation
% 
% figure
% clc
% [vMat, fMat] = spheretri(30000);
% patch('Vertices',vMat,'Faces',fMat,'FaceColor','g','EdgeColor','k');
% 
% % quick check on points for first face..
% face = vMat(fMat(1,:), :);
% [flon, flat, fr] = cart2sph(face(:,1), face(:,2), face(:,3) );
% flat = rad2deg(flat);
% flon = rad2deg(flon);
% 
% s1 = m_idist(flon(1), flat(1), flon(2), flat(2));
% s2 = m_idist(flon(1), flat(1), flon(3), flat(3));
% s3 = m_idist(flon(2), flat(2), flon(3), flat(3));
% 
% % distance in degrees between vertices on particular face..
% fd = [km2deg(s1./1e3), km2deg(s2./1e3), km2deg(s3./1e3)];
% 
% 
% [llon, llat, rr] = cart2sph(vMat(:,1), vMat(:,2), vMat(:,3) );
% llat = rad2deg(llat); llon = rad2deg(llon);
% 
% 
% 
% %% 6. select faces located exclusively in africa
% 
% %
% clc
% nFaces = length(fMat);
% faceInAfrica = false(nFaces, 1);  % about 82k points for 1deg resl.
% 
% for iFace = 1 : nFaces
%     face = vMat(fMat(iFace,:), :);
%     [flon, flat, ~] = cart2sph(face(:,1), face(:,2), face(:,3) );
%     flat = rad2deg(flat); mflat = mean(flat);
%     flon = rad2deg(flon); mflon = mean(flon);
%     
%     INaf = inpolygon(flon,flat,aflon,aflat);
%     INma = inpolygon(flon,flat,Malon,Malat);
%     numIN = floor(sum(INaf)) + floor(sum(INma));
%     
%     % if all 3 vertices are inside africa then face is in!
% %     if numIN == 3
% %%%%% Siyu: now try if one one vertix is sufficient
%     if numIN >= 1
%         faceInAfrica(iFace) = true;
%         disp([flat, flon]); %CHECKED all points are in Hawaii
%         %disp(iFace)
%     end
%     
% end
% 
% afMat = fMat(faceInAfrica,:);
% 
% % unique vertices
% %afVI = unique(afMat, 'rows'); % unique index of all faces in africa ...
% afVI = find(faceInAfrica == 1);
% % afVMat = vMat(afMat,:);
% 
% % [allon, allat, ~] = cart2sph(afVMat(:,1), afVMat(:,2), afVMat(:,3) );
% % allat = rad2deg(allat); allon = rad2deg(allon);
% 
% %% 7 KERNEL: on what face do points lie ?
% clc
% 
% % set the distance threshold for the three plots
% Degshort = 10;
% Degmid = 20;
% Deglong = 30;
% 
% nFaces = length(afMat);
% cntInFace = zeros(nFaces, 1);  % about 82k points for 1deg resl.
% 
% faceHit = 0;
% faceInd = [];   % index of face with 
% %faceCnt = [];
% 
% faceStat = zeros(nFaces, 3); % center coordinate of all faces, with count of points on face
% faceStatL = zeros(nFaces, 3); % number of faces with path greater than the shord, mid, long threshold
% 
% [totPath, ~] = size(ploc) ;
% 
% for iK = 1:totPath
%     ln = [ploc(iK,2); ploc(iK,4)];
%     lt = [ploc(iK,3); ploc(iK,5)];
%     
%     [s, az12, az21] = m_idist(ln(1), lt(1), ln(2), lt(2));
%     sdel = km2deg(s./1e3);
%     
%     % points on path at 0.25 degree resolution
%     nlegs = floor(sdel/delP);
%     [ltp, lnp] = gcwaypts(lt(1), ln(1), lt(2), ln(2), nlegs);
%     
%     siK = num2str(iK); snK = num2str(totPath);
%     nFace4Path = 0;
% 
%     for iFace = 1 : length(afVI)
%         face = vMat(afMat(iFace,:), :);
%         [flon, flat, fr] = cart2sph(face(:,1), face(:,2), face(:,3) );
%         flat = rad2deg(flat); mflat = mean(flat);
%         flon = rad2deg(flon); mflon = mean(flon);
%         
%         
%         %  see if path is in face ..
%         IN = inpolygon(lnp,ltp,flon,flat);
%         numIN = floor(sum(IN));
%         % disp(numIN);   % seems that none of the path is in any face
%         %[iK iFace]
%         if numIN >= 1
%             %faceStat(iFace, 3) = [mflon, mflat, numIN];
%             siFace = num2str(iFace);
%             faceStat(iFace, 3) = faceStat(iFace, 3) + 1;
% 
%             % Dr. O's old version of the plots     
% %             if sdel <= maxDeg
% %                 faceStatL(iFace,1) = faceStatL(iFace,1) + 1;
% %             else
% %                 faceStatL(iFace,2) = faceStatL(iFace,2) + 1;
% %             end
% 
%             if sdel <= Degshort
%                 faceStatL(iFace,1) = faceStatL(iFace,1) + 1;
%             elseif sdel <= Degmid
%                 faceStatL(iFace,2) = faceStatL(iFace,2) + 1;
%             elseif sdel <= Deglong
%                 faceStatL(iFace,3) = faceStatL(iFace,3) + 1;
%             end
% 
%             nFace4Path = nFace4Path + 1;
%         end   
%         
%     end
%     
%     sFace4Path = num2str(nFace4Path);
%     snFaces = num2str(nFaces);
%     
%     disp(['completed, path ' siK ' of ' snK ' matched ' sFace4Path ...
%         ' of ' snFaces]);
% end


%% visualize ray coverage ...
% load the raypath data
load('/Users/sxue3/Documents/BayMap_Figures/Data/ADAMA_raypath.mat')

figure
clf

latlim = [-39 40];
lonlim = [-26 60];

for iFace = 1 : length(afVI)
    face = vMat(afMat(iFace,:), :);
    [flon, flat, fr] = cart2sph(face(:,1), face(:,2), face(:,3) );
    flat = rad2deg(flat); mflat = mean(flat);
    flon = rad2deg(flon); mflon = mean(flon);
    
    faceStat(iFace, 1:2) = [mflon, mflat];

end

% visualize count ... 
sz = .25;
lond = -26:sz:60;
latd = -39:sz:40;

[dlon, dlat] = meshgrid(lond, latd);

allcount = faceStat(:,3);
allcount(allcount==0)=NaN;

F = scatteredInterpolant(faceStat(:,1), faceStat(:,2), allcount); 
F.ExtrapolationMethod = 'none';

fval = F(dlon, dlat) ;
fval(isnan(fval))=0;

% Create a mask to seperate region within Africa
afrmask = zeros(317, 345);
mapsize = size(afrmask);
for xi = 1:mapsize(2)
    loni = dlon(1,xi);
    for yi = 1:mapsize(1)
        lati = dlat(yi,1);

        if isinterior(coastline, loni, lati) == 1 || isinterior(madgline, loni, lati) == 1
            afrmask(yi, xi) = 1;
        end
    end
end

afrmask = logical(afrmask - 1);  % points not in Africa == 1

% apply the African mask
fval(afrmask) = nan;

     
% -- plot face count

m_proj('miller', 'lat',latlim, 'lon', lonlim );
hold on
m_grid('box','fancy','tickdir','in');
hold on
m_pcolor(dlon, dlat, log10(fval))
cb = colorbar('eastoutside', 'fontsize', 15);
caxis([min(log10(fval),[],'all'), max(log10(fval),[],'all')])
colormap(flipud(jet))
ylabel(cb, 'log_{10}[# Rays]', 'fontsize', 20);
m_line(aflon, aflat, 'color','k','linewi',2);
m_line(Malon, Malat, 'color','k','linewi',2);

% -- cratons
load('/Users/sxue3/Documents/BayMap_Figures/Data/AfricanCratons.mat');
% plot cratons
m_line(AC.Congo(:,2), AC.Congo(:,1), 'color',[0.3, 0.3, 0.3],'linewi',2);
m_line(AC.Kala(:,2), AC.Kala(:,1), 'color',[0.3, 0.3, 0.3],'linewi',2);
m_line(AC.Sahara(:,2), AC.Sahara(:,1), 'color',[0.3, 0.3, 0.3],'linewi',2);
m_line(AC.Tanza(:,2), AC.Tanza(:,1), 'color',[0.3, 0.3, 0.3],'linewi',2);
m_line(AC.West(:,2), AC.West(:,1), 'color',[0.3, 0.3, 0.3],'linewi',2);
m_line(AC.Zim(:,2), AC.Zim(:,1), 'color',[0.3, 0.3, 0.3],'linewi',2);

% Plot the down sampled points
load('/Users/sxue3/Documents/BayMap_Figures/Data/AfrDownSampleDots.mat');

for i = 1:length(dslat)

     m_line(dslon(i), dslat(i), 'marker','^','color','k','linewi',1,...
        'linest','none','markersize',1,'markerfacecolor','w');
end


% set the size of the plot
x0=10;
y0=10;
width=800;
height=800;
set(gcf,'position',[x0,y0,width,height]);

%%%% Plot the litho Raypath coverage
ax2 = axes('position',[0.134, 0.15, 0.3, 0.35]);
ax2.YAxis.Visible = 'off'; % remove y-axis
ax2.XAxis.Visible = 'off'; % remove x-axis

load('/Users/sxue3/Documents/BayMap_Figures/Data/Litho_raypath.mat')
m_proj('lambert','lon',[-26 60],'lat',[-39 40]);  % African boundaries
m_pcolor(LG,LT,InterpDat);
hold on
m_grid('linestyle','none','xtick',[],'ytick',[],'xticklabel',[],'yticklabel',[])
m_line(aflon, aflat, 'color','k','linewi',2);
m_line(Malon, Malat, 'color','k','linewi',2);

% plot cratons
m_line(AC.Congo(:,2), AC.Congo(:,1), 'color',[0.8, 0.8, 0.8],'linewi',1.5);
m_line(AC.Kala(:,2), AC.Kala(:,1), 'color',[0.8, 0.8, 0.8],'linewi',1.5);
m_line(AC.Sahara(:,2), AC.Sahara(:,1), 'color',[0.8, 0.8, 0.8],'linewi',1.5);
m_line(AC.Tanza(:,2), AC.Tanza(:,1), 'color',[0.8, 0.8, 0.8],'linewi',1.5);
m_line(AC.West(:,2), AC.West(:,1), 'color',[0.8, 0.8, 0.8],'linewi',1.5);
m_line(AC.Zim(:,2), AC.Zim(:,1), 'color',[0.8, 0.8, 0.8],'linewi',1.5);


colormap(flipud(jet))
caxis([min(log10(fval),[],'all'), max(log10(fval),[],'all')])

% set(gca,'xticklabel',[])
ax2.YAxis.Visible = 'off'; % remove y-axis
ax2.XAxis.Visible = 'off'; % remove x-axis

% label the maps
annotation('textbox',[.7 .63 .4 .2], 'String','ADAMA','EdgeColor','none', 'FontSize', 20)
annotation('textbox',[.15 .05 .4 .2], 'String','Litho1.0','EdgeColor','none', 'FontSize', 20)

% save the plot in PDF format
fig = gcf;
saveFig('fig2_RayPathDensity.pdf', './fig/', 1, fig);

