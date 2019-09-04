function DisplayData(obs)

% Script to Display Data

%% Pre-flight

% figure defaults
DGdisplaydefaults;

if ~exist('obs','var')
    clear, clc, close all
    obs = 'LM';
end

cols = jet(16);

%% Load Data

data_folder = fullfile('C:','Users','cege-user','Dropbox','UCL','Data','LargeSphere','Experimental Data');
if strcmp(obs,'LM')
    rootdir = fullfile(data_folder,'\2012 Sep LM\Results - Sep 2012');    
    dfile = 'C:\Users\cege-user\Dropbox\UCL\Data\LargeSphere\Hardware Data\Calibration Data\LCD display measurement.mat';
elseif strcmp(obs,'TR')
    rootdir = fullfile(data_folder,'\2013 Apr TR');
    dfile = 'C:\Users\cege-user\Dropbox\UCL\Data\LargeSphere\Hardware Data\Calibration Data\Large LCD display measurement.mat'; % no characterisation data matching the XYZ listed in the AIC paper seems to be available. This one is close, and could potentially have been made at around the right time.
elseif strcmp(obs,'DG')
    rootdir = fullfile(data_folder,'\2016 Oct DG');
    dfile = 'C:\Users\cege-user\Dropbox\UCL\Data\LargeSphere\Hardware Data\Calibration Data\Large LCD display measurement - Oct 2016.mat';
elseif strcmp(obs,'baseline')
    rootdir = fullfile(data_folder,'\BaselineData');
    dfile = 'C:\Users\cege-user\Dropbox\UCL\Data\LargeSphere\Hardware Data\Calibration Data\Large LCD display measurement - Oct 2016.mat'; % This is an arbitrary pick. Really I should do it with each of the dfiles listed above
else
    error('No/incorrect observer entered')
end

cd(rootdir)
files = dir('*nm*.mat');
if strcmp(obs,'baseline')
    files = dir('*.mat');
end
    
for j=1:length(files)    
    load(fullfile(rootdir,files(j).name));  % load experimental results
    files(j).dataLAB = LABmatch;
    %LAB(:,:,:,j) = LABmatch;
    files(j).dataRGB = RGBmatch;
    files(j).Tmatch  = Tmatch;    
end

LN = size(files(1).dataLAB,2);      % number of lightness levels per repeat
N = size(files(1).dataLAB,3);       % number of repetitions over time


%% Remove repeats

if strcmp(obs,'LM')
    files(4) = [];
    files(5) = [];
    files(6) = [];
    files(6) = []; % repeating numbers because each one is removed one-by-one and the rest drop back by one to fill the gap
    files(6) = [];
    %LAB(:,:,:,[4,6,8:10]) = []; %different number because we're removing them simultaneously
end

if strcmp(obs,'TR')
    files(4) = [];
    files(6) = [];
    %LAB(:,:,:,[4,7]) = [];
end

if strcmp(obs,'DG')
    files(6) = [];
    files(6) = [];
    files(6) = [];
    %LAB(:,:,:,6:8) = [];
end

%% Creates calibrated LAB values

%load calibration file
load(dfile,'sval','XYZ')

%interpolate recorded values (sval) to required vals (0:1:255)
XYZinterp=zeros(3,256,4);
% figure, hold on
for i=1:3
    for j=1:4
        XYZinterp(i,:,j) = interp1(sval, XYZ(i,:,j), 0:255, 'spline');
%         plot(sval, XYZ(i,:,j),'o')        % check interpolation 
%         plot(0:255, XYZinterp(i,:,j))
    end
end

% Calcaulate XYZ for white point of display
%   This method gives slightly different results to the previous method
%   (where cie1931 was loaded, and fresh tristimulus were calculated from
%   the recorded spectra, but this method is much neater and in-ilne with
%   the use of the PR650 XYZ values elsewhere).

XYZw = XYZ(:,end,4)/XYZ(2,end,4)*100;

for trial = 1:length(files)   

    % Thresholding:
    %   Original RGB values included out of gamut (sRGB)
    %   selections, resulting in above 1 and below 0 values. These would
    %   actually have only presented values at 0/1 and so here they are
    %   corrected to represent what would actually have been presented    
    
    files(trial).dataRGBgamflag = files(trial).dataRGB > 1 | files(trial).dataRGB < 0; %out of gamut flag
    
    files(trial).dataRGBgamcor  = files(trial).dataRGB; %duplicate
    files(trial).dataRGBgamcor(files(trial).dataRGB < 0) = 0;
    files(trial).dataRGBgamcor(files(trial).dataRGB > 1) = 1;
    
    % Quantization
    files(trial).dataRGBgamcor = uint8(files(trial).dataRGBgamcor*255);
    
    files(trial).dataXYZcal    = zeros(3,LN,N);
    files(trial).dataxycal     = zeros(2,LN,N);
    files(trial).dataLABcal    = zeros(3,LN,N);
    
    for j = 1:LN
        for k = 1:N
            files(trial).dataXYZcal(:,j,k) = ...
                (XYZinterp(:,files(trial).dataRGBgamcor(1,j,k)+1,1)...
                +XYZinterp(:,files(trial).dataRGBgamcor(2,j,k)+1,2)...
                +XYZinterp(:,files(trial).dataRGBgamcor(3,j,k)+1,3));
            
            files(trial).dataxycal(1,j,k)=...
                files(trial).dataXYZcal(1,j,k)/sum(files(trial).dataXYZcal(:,j,k));
            files(trial).dataxycal(2,j,k)=...
                files(trial).dataXYZcal(2,j,k)/sum(files(trial).dataXYZcal(:,j,k));
            
            files(trial).dataLABcal(:,j,k)=...
                XYZToLab(files(trial).dataXYZcal(:,j,k),XYZw);
        end
    end
end

% % Plot display white points
figure, hold on
drawChromaticity
scatter(XYZw(1)/sum(XYZw),XYZw(2)/sum(XYZw),'k')

%%

for i = 1:length(files)
    LAB(:,:,:,i) = files(i).dataLABcal;
end


%% Basic scatter of everything

figure, hold on
scatter3(LAB(2,:),LAB(3,:),LAB(1,:),'.') %how would I add colour to this?

xlabel('a*')
ylabel('b*')
zlabel('L*')

%% Grouped by repeat and wavelength

cla

for j=1:size(LAB,4)
    for i=1:size(LAB,3)
        plot3(LAB(2,:,i,j),LAB(3,:,i,j),LAB(1,:,i,j),'Color',cols(j,:))
    end
end

%% Grouped by wavelength, averaged over repeats

figure('Position',[100 100 1000 800])
hold on
subplot(2,2,1)
hold on

xlabel('a*')
ylabel('b*')
zlabel('L*')

LABm = squeeze(median(LAB,3));

for j=1:size(LAB,4)   
    p3(j) = plot3(LABm(2,:,j),LABm(3,:,j),LABm(1,:,j),'o-',...
        'Color',cols(j,:),'DisplayName',files(j).name(1:5));    
end

view(2)
legend('Location','westoutside')

s2 = subplot(2,2,3);
xlabel('a*')
ylabel('b*')
zlabel('L*')
copyobj(p3,s2);
view(0,0)

s3 = subplot(2,2,4);
xlabel('a*')
ylabel('b*')
zlabel('L*')
copyobj(p3,s3);
view(90,0)


%% Select Lstar
subplot(2,2,2)
hold on

%ltns_all = 85:-5:10;
%ltns = [16,9,1]; %lightnesses
ltns = [14,6]; %lightnesses
ltnscols = [0,0,0; 0.5,0.5,0.5;0.7,0.7,0.7];
linestyle = {'-','--',':'};

for lI = 1:length(ltns) %lightness Index
     p(lI) = plot3(squeeze(LABm(2,ltns(lI),:)),squeeze(LABm(3,ltns(lI),:)),squeeze(LABm(1,ltns(lI),:)),...
         'Color',ltnscols(lI,:),'LineStyle',linestyle{lI});
    for j=1:size(LAB,4)
        scatter3(LABm(2,ltns(lI),j),LABm(3,ltns(lI),j),LABm(1,ltns(lI),j),...
            [],cols(j,:),'filled')
    end
end

%legend(p,split(num2str(ltns_all(ltns)))) %programmatically
%legend(p,{'10 L*','45 L*','85 L*'},'Location','best')
legend(p,{'20 L*','60 L*'},'Location','best')
view(2)
xlabel('a*')
ylabel('b*')
zlabel('L*')

%save2pdf([data_folder(1:end-17),'Data Analysis\figs\',obs,'dataOverview'])

%%

% for k = 1:16
%     for i=1:16
%         for j=1:16
%             [H(i,j,k), pValue(i,j,k)] = kstest_2s_2d(squeeze(LAB(1:2,k,:,i))', squeeze(LAB(1:2,k,:,j))');
%         end
%     end
% end
% 
% figure, imagesc(mean(pValue,3))
% colorbar

end