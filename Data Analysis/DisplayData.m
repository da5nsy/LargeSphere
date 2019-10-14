function DisplayData(obs)

% Script to Display Data

%% Pre-flight

% figure defaults
DGdisplaydefaults;

if ~exist('obs','var')
    clear, clc, close all
    obs = 'LM';
end

cs = 'LAB'; %colour space, for plotting
cols = jet(16);
%load cols cols
bcol = lines(9);

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
XYZinterp = zeros(3,256,4);
% figure, hold on
for i = 1:3
    for j = 1:4
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
            
            files(trial).dataxycal(1,j,k) = ...
                files(trial).dataXYZcal(1,j,k)/sum(files(trial).dataXYZcal(:,j,k));
            files(trial).dataxycal(2,j,k) = ...
                files(trial).dataXYZcal(2,j,k)/sum(files(trial).dataXYZcal(:,j,k));
            files(trial).dataxycal(3,j,k) = files(trial).dataXYZcal(2,j,k);
            
            files(trial).dataLABcal(:,j,k) = ...
                XYZToLab(files(trial).dataXYZcal(:,j,k),XYZw);
        end
    end
end

% % Plot display white points
figure, hold on
DrawChromaticity
scatter(XYZw(1)/sum(XYZw),XYZw(2)/sum(XYZw),'k')

%% Grab LABcal values out for easier plotting

for i = 1:length(files)
    LAB(:,:,:,i) = files(i).dataLABcal;
    xyY(:,:,:,i) = files(i).dataxycal;
end


%% Basic scatter of everything

if strcmp(cs,'LAB')
    figure, hold on
    scatter3(LAB(2,:),LAB(3,:),LAB(1,:),'.') %how would I add colour to this?
    
    xlabel('a*')
    ylabel('b*')
    zlabel('L*')
    
elseif strcmp(cs,'xyY')
    figure, hold on
    DrawChromaticity
    scatter3(xyY(1,:),xyY(2,:),xyY(3,:),'.') %how would I add colour to this?
    
    xlabel('x')
    ylabel('y')
    zlabel('Y')
end

%% Grouped by repeat and wavelength

cla

if strcmp(cs,'LAB')
    for j = 1:size(LAB,4)
        for i=1:size(LAB,3)
            if ~strcmp(obs,'baseline')
                plot3(LAB(2,:,i,j),LAB(3,:,i,j),LAB(1,:,i,j),'Color',cols(j,:))
            else
                plot3(LAB(2,:,i,j),LAB(3,:,i,j),LAB(1,:,i,j),'--','Color',bcol(ceil(j/5),:))
            end
        end
    end
elseif strcmp(cs,'xyY')
    DrawChromaticity
    daspect([1,1,150])
    for j = 1:size(xyY,4)
        for i = 1:size(xyY,3)
            if ~strcmp(obs,'baseline')
                plot3(xyY(1,:,i,j),xyY(2,:,i,j),xyY(3,:,i,j),'Color',cols(j,:))
            else
                plot3(xyY(1,:,i,j),xyY(2,:,i,j),xyY(3,:,i,j),'--','Color',bcol(ceil(j/5),:))
            end
        end
    end
end

%% Grouped by wavelength, averaged over repeats

figure('Position',[100 100 1000 800])
hold on
if strcmp(obs,'baseline')
    subplot(2,2,1:2)
    hold on
else
    subplot(2,2,1)
    hold on
end

if strcmp(cs,'xyY')
    xlabel('x')
    ylabel('y')
    zlabel('Y')
    
    xyYm = squeeze(median(xyY,3));
    
    for j=1:size(xyY,4)
        p3(j) = plot3(xyYm(1,:,j),xyYm(2,:,j),xyYm(3,:,j),'o-',...
            'Color',cols(j,:),'DisplayName',files(j).name(1:5));
        daspect([1,1,150])
    end
    
    view(2)
    legend('Location','westoutside')
    
    s2 = subplot(2,2,3);
    xlabel('x')
    ylabel('y')
    zlabel('Y')
    copyobj(p3,s2);
    view(0,0)
    
    s3 = subplot(2,2,4);
    xlabel('x')
    ylabel('y')
    zlabel('Y')
    copyobj(p3,s3);
    view(90,0)
    
    subplot(2,2,2)
    hold on
    
    %ltns_all = 85:-5:10;
    %ltns = [16,9,1]; %lightnesses
    ltns = [14,6]; %lightnesses
    ltnscols = [0.5,0.5,0.5;0.7,0.7,0.7;0,0,0];
    linestyle = {'--',':','-'};
    
    for lI = 1:length(ltns) %lightness Index
        p(lI) = plot3(squeeze(xyYm(1,ltns(lI),:)),squeeze(xyYm(2,ltns(lI),:)),squeeze(xyYm(3,ltns(lI),:)),...
            'Color',ltnscols(lI,:),'LineStyle',linestyle{lI});
        for j=1:size(LAB,4)
            scatter3(xyYm(1,ltns(lI),j),xyYm(2,ltns(lI),j),xyYm(3,ltns(lI),j),...
                [],cols(j,:),'filled')
        end
    end
    
    % averaging option
    %p = plot3(squeeze(mean(xyYm(1,6:14,:),2)),squeeze(mean(xyYm(2,6:14,:),2)),squeeze(mean(xyYm(3,6:14,:),2)),'k');
    
    
    %legend(p,split(num2str(ltns_all(ltns)))) %programmatically
    %legend(p,{'10 L*','45 L*','85 L*'},'Location','best')
    %legend(p,{'20 L*','60 L*'},'Location','best')
    view(2)
    xlabel('x')
    ylabel('y')
    zlabel('Y')
    
elseif strcmp(cs,'LAB')
    
    xlabel('a*')
    ylabel('b*')
    zlabel('L*')
    
    LABm = squeeze(median(LAB,3));
    
    t = 3:5:length(files); %[3,8,13,18,23];
    if strcmp(obs,'baseline')
        for j=1:length(t)%1:size(LAB,4)
            p3(j) = plot3(LABm(2,:,t(j)),LABm(3,:,t(j)),LABm(1,:,t(j)),'o-',...
                'Color',bcol(ceil(t(j)/5),:),'DisplayName',files(t(j)).name(1:end-10));
        end
        legend('Location','westoutside')
    else
        for j=1:size(LAB,4)
            p3(j) = plot3(LABm(2,:,j),LABm(3,:,j),LABm(1,:,j),'o-',...
                'Color',cols(j,:),'DisplayName',files(j).name(1:5));
        end
        legend('Location','westoutside')
    end
    
    if strcmp(obs,'baseline')
        axis equal
        xlim([-100 100])
        ylim([-100 100])
        zlim([0 100])
        cleanTicks
    end
    view(2)
    
    s2 = subplot(2,2,3);
    xlabel('a*')
    ylabel('b*')
    zlabel('L*')
    copyobj(p3,s2);
    view(0,0)
    if strcmp(obs,'baseline')
        axis equal
        xlim([-100 100])
        ylim([-100 100])
        zlim([0 100])
        cleanTicks
    end
    
    s3 = subplot(2,2,4);
    xlabel('a*')
    ylabel('b*')
    zlabel('L*')
    copyobj(p3,s3);
    view(90,0)
    if strcmp(obs,'baseline')
        axis equal
        xlim([-100 100])
        ylim([-100 100])
        zlim([0 100])
        cleanTicks
    end
    
    if ~strcmp(obs,'baseline')
        subplot(2,2,2)
        hold on
    end
    
    %ltns_all = 85:-5:10;
    %ltns = [16,9,1]; %lightnesses
    ltns = [14,6]; %lightnesses
    ltnscols = [0.5,0.5,0.5;0.7,0.7,0.7;0,0,0];
    linestyle = {'--',':','-'};
    
    for lI = 1:length(ltns) %lightness Index
        if strcmp(obs,'baseline')
            %             lI = 2;
            %             for j = [t,t(1)] %repeat initial (central)
            %                 for i = [-2,-1,1,2]
            %                     p(lI) = plot3([LABm(2,ltns(lI),j),LABm(2,ltns(lI),j+i)],[LABm(3,ltns(lI),j),LABm(3,ltns(lI),j+i)],[LABm(1,ltns(lI),j),LABm(1,ltns(lI),j+i)],...
            %                         'Color',ltnscols(lI,:),'LineStyle',linestyle{lI});
            %                 end
            %             end
            %             legend(p,{'20 L*','60 L*'},'Location','best')
            %             xlim([-100 100])
            %             ylim([-100 100])
            %             zlim([0 100])
            %             axis equal
            %             cleanTicks
            %             for j=1:size(LAB,4)
            %                 scatter3(LABm(2,ltns(lI),j),LABm(3,ltns(lI),j),LABm(1,ltns(lI),j),...
            %                     [],bcol(ceil(j/5),:),'filled')
            %             end
            %             view(2)
            %             xlabel('a*')
            %             ylabel('b*')
            %             zlabel('L*')
        else
            p(lI) = plot3(squeeze(LABm(2,ltns(lI),:)),squeeze(LABm(3,ltns(lI),:)),squeeze(LABm(1,ltns(lI),:)),...
                'Color',ltnscols(lI,:),'LineStyle',linestyle{lI});
            %legend(p,split(num2str(ltns_all(ltns)))) %programmatically
            %legend(p,{'10 L*','45 L*','85 L*'},'Location','best')
            legend(p,{'20 L*','60 L*'},'Location','best','AutoUpdate','off')
            for j=1:size(LAB,4)
                scatter3(LABm(2,ltns(lI),j),LABm(3,ltns(lI),j),LABm(1,ltns(lI),j),...
                    [],cols(j,:),'filled')
            end
            view(2)
            xlabel('a*')
            ylabel('b*')
            zlabel('L*')
        end
    end
end

%save2pdf([data_folder(1:end-17),'Data Analysis\figs\',obs,'dataOverview'])

%%

figure('Position',[100 100 500 800])
hold on

for i = 1:3    
    subplot(3,1,i)
    imagesc(squeeze(median(LAB(i,6:end,:,:),2)))
    %imagesc(squeeze(LAB(i,6,:,:)))
    axis image
    
    yticks(1:10)
    
    ylabel('Repeat #')
    colormap('gray')
    colorbar
    
    if i == 3
        xticks(1:16)
    xticklabels(400:20:700)
    xtickangle(90)
    xlabel('Wavelength (nm)')
    else
        xticks('')
    end
    set(gca,'YDir','normal')
end

%save2pdf([data_folder(1:end-17),'Data Analysis\figs\',obs,'dataOverTime'])

%%
if strcmp(obs,'LM')
    load adaptingFieldCIELabLM.mat Lab
elseif strcmp(obs,'TR')
    load adaptingFieldCIELabTR.mat Lab
else
    error('I haven''t made a set of surround CIELAB values for that observer')
end
adap = 400:20:700;

figure, hold on

plot3(Lab(2,:),Lab(3,:),Lab(1,:),'k')
for i = 1:16
    scatter3(Lab(2,i),Lab(3,i),Lab(1,i),[],cols(i,:),'filled','DisplayName',num2str(adap(i)))
end
grid off
xlabel('a*')
ylabel('b*')
zlabel('L*')

% for j=1:size(LAB,4)
%     p3(j) = plot3(LABm(2,:,j),LABm(3,:,j),LABm(1,:,j),'o-',...
%         'Color',cols(j,:),'DisplayName',files(j).name(1:5));
% end
%legend('Location','westoutside')

for lI = 1:length(ltns) %lightness Index
    p(lI) = plot3(squeeze(LABm(2,ltns(lI),:)),squeeze(LABm(3,ltns(lI),:)),squeeze(LABm(1,ltns(lI),:)),...
        'Color',ltnscols(lI,:),'LineStyle',linestyle{lI});
    for j=1:size(LAB,4)
        scatter3(LABm(2,ltns(lI),j),LABm(3,ltns(lI),j),LABm(1,ltns(lI),j),...
            [],cols(j,:),'filled')
    end
end

 save2pdf([data_folder(1:end-17),'Data Analysis\figs\',obs,'compareWithSurround'])


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