function DisplayData(obs)

% Script to Display Data

%% Pre-flight

% figure defaults
set(groot,'defaultfigureposition',[100 100 500 400]);
set(groot,'defaultLineLineWidth',2);
set(groot,'defaultAxesFontName', 'Courier');
set(groot,'defaultAxesFontSize',12);
set(groot,'defaultFigureRenderer', 'painters') %renders pdfs as vectors
set(groot,'defaultfigurecolor','white')

if ~exist('obs','var')
    clear, clc, close all
    obs = 'baseline';
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
    disp('No/incorrect observer entered')
end

cd(rootdir)
files = dir('*nm*.mat');
if strcmp(obs,'baseline')
    files = dir('*.mat');
end
    

%N = 10;                             % number of repetitions over time
%LN = 16;                            % number of lightness levels per repeat

for j=1:length(files)    
    load(fullfile(rootdir,files(j).name));  % load experimental results
    files(j).dataLAB = LABmatch;
    LAB(:,:,:,j) = LABmatch;
    files(j).dataRGB = RGBmatch;
    files(j).Tmatch  = Tmatch;    
end

%% Remove repeats

if strcmp(obs,'LM')
    files(4) = [];
    files(5) = [];
    files(6) = [];
    files(6) = []; % repeating numbers because each one is removed one-by-one and the rest drop back by one to fill the gap
    files(6) = [];
    LAB(:,:,:,[4,6,8:10]) = []; %different number because we're removing them simultaneously
end

if strcmp(obs,'TR')
    files(4) = [];
    files(6) = [];
    LAB(:,:,:,[4,7]) = [];
end

if strcmp(obs,'DG')
    files(6) = [];
    files(6) = [];
    files(6) = [];
    LAB(:,:,:,[6:8]) = [];
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

LABm = squeeze(mean(LAB,3));

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