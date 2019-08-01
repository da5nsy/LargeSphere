% I've grabbed some code from Analyse results by cone excitation vs wavelength.m
% I want to get it to run with demo data, to show what the case would be if
% it was a simple Von Kries scaling.

%%

clear, clc, close all

plt = 1;

% Load observer
load T_cones_ss10.mat
load T_rods.mat
load T_melanopsin.mat
T_obs = [T_cones_ss10;...
    SplineCmf(S_rods,T_rods,S_cones_ss10);...
    SplineCmf(S_melanopsin,T_melanopsin,S_cones_ss10)];
S_obs = S_cones_ss10;
clear T_cones_ss10 S_cones_ss10 T_rods S_rods T_melanopsin S_melanopsin %cleanup
%figure, plot(SToWls(S_obs),T_obs')

%% Grabbing bits to make an empty data container

% Load spectrums of peripheral stimuli
load('C:\Users\cege-user\Dropbox\UCL\Data\LargeSphere\Hardware Data\Filter spectra\Illumination in sphere.mat','spectra')
if plt
    figure, plot(SToWls([380,4,101]),spectra(:,2:17))
end
spectra = SplineSpd([380,4,101],spectra(:,2:17),S_obs);

wmin = 400;                         % minimum wavelength filter
wmax = 700;                         % maximum wavelength filter
wrange = wmin:20:wmax;              % range of wavelengths (20nm intervals)
N = 10;                             % number of repetitions over time
LN = 16;                            % number of lightness levels per repeat
WN = length(wrange);                % number of wavelength samples
Lval = 85:-5:10;                    % lightness values

%TNM = floor(min(TI(4,:)/60));      % find length of shortest session (minutes)
TNM = 72;                           % (DG note: This is what the real 2013 data ends up as)

LMSR_sim = zeros(size(T_obs,1),LN,TNM,WN,'double'); % Empty data container (Should be 4,16,72,16)

% SO, what do these number mean?

% 4 - Long, Medium, Short, Rod, Mel
% 16 - number of lightness levels (LN)
% 72 - minutes of session
% 16 - number of wavelength samples

%% Generate simulated data

norm = 1;

for ii = 1:size(LMSR_sim,1)    % each sensor
    for jj = 1:LN              % (let's just do one lightness level for now)
        for kk = 1:TNM         % (let's assume the same across time for now)
            for mm = 1:WN      % for each wavelength sample
                
                light = spectra(:,mm);
                sensor = T_obs(ii,:);
                stimulation = sensor * light;
                nsfwp = norm * stimulation; %needed stimulation for white point
                LMSR_sim(ii,jj,kk,mm) = nsfwp;
            end
        end
    end
end

%% Plot 3D surfaces of LMSR vs wavelength and time for constant lightness

lness = 6:11;              % lightness (can be single or multiple)
LMSI_sim=squeeze(mean(LMSR_sim(:,lness,:,:),2));
labels = {'L','M','S','R','I'};

if plt
    figure('units','normalized','outerposition',[0 0 1 1])
    for i=1:size(T_obs,1)
        subplot(1,size(T_obs,1),i)
        imagesc(squeeze(LMSI_sim(i,:,:)))
        set(gca,'YDir','normal')
        %colorbar
        
        xticks(1:16)
        xticklabels(wrange)
        xlabel({'Wavelength of adapting field (nm)',labels{i}});
        if i == 1
            ylabel('Time (min)');
        end
        colormap gray
    end
end

%% Load a real dataset
load 'C:\Users\cege-user\Dropbox\UCL\Data\LargeSphere\Data Analysis\Tania LMSI for testing VonKriesTest'
LMSI_real = LMSI; clear LMSI

if plt
    figure('units','normalized','outerposition',[0 0 1 1])
    for i=1:3
        subplot(1,3,i)
        imagesc(squeeze(LMSI_real(i,:,:)))
        set(gca,'YDir','normal')
        %colorbar
        
        xticks(1:16)
        xticklabels(wrange)
        xlabel({'Wavelength of adapting field (nm)',labels{i}});
        if i == 1
            ylabel('Time (min)');
        end
        colormap gray
    end
end

%% Try random combos of the cones to see if we can emulate the real data

rng(1)

for q = 1:100 % Run it X number of times to see whether difference between different values of nIn is just noise or something meaningful
    
    nIn = 5;            %number of inputs (L,M,S,R,I in that order)
    nOut = 3;
    runs = 1000;
    SFM = 50;           %Scaling factor max
    
    for i=1:nIn % This creates duplicate data and no longer needs to be a cell (previously I thresholded the data that I stored based on correlation which left different columns with different lengths, but now we keep everything.
        CWstore{i} = zeros(nIn,runs);  %Cone weight store
    end
    
    for i=1:nOut
        Cstore{i} = 0;              %Correlation store
    end
    
    for j=1:runs
        ConeWeights = rand(nIn,1)*SFM-(SFM/2);
        randomEffect = zeros([TNM,WN]);
        for i=1:nIn
            randomEffect = randomEffect + squeeze((LMSI_sim(i,:,:)*ConeWeights(i,1)));
        end
        
        %     figure
        %     imagesc(randomEffect)
        %     set(gca,'YDir','normal')
        %     xticks(1:16)
        %     xticklabels(wrange)
        %     colormap gray
        
        randomEffect_slim = mean(randomEffect,1);
        
        % Compare random effect with our actual data
        LMSI_real_slim = squeeze(mean(LMSI_real,2));
        for k=1:nIn
            CWstore{k} = cat(2,CWstore{k},ConeWeights);
        end
        for k=1:nOut % Just trying to predict L,M,S (for now)
            
            Cstore{k} = cat(2,Cstore{k},...
                min(min(corrcoef(randomEffect_slim,LMSI_real_slim(k,:)'))));
        end
        if mod(j,10000) == 0 % to read out progress when running large runs
            disp(j)
        end
    end
    
    %figure,
    for i=1:nOut
        %subplot(1,nOut,i)
        %histogram(abs(Cstore{1,i}))
        %xlim([0 1])
        mc(i,q) = max(abs(Cstore{1,i}));
    end
end

disp(mean(mc,2)')
disp(max(mc'))

%% Visualise best performance

figure,
for k=1:nOut
    subplot(2,nOut,k)
    [~,j] = max(abs(Cstore{1,k}));
    randomEffect = zeros([TNM,WN]);
    for i=1:nIn
        randomEffect = randomEffect + squeeze((LMSI_sim(i,:,:)*CWstore{1,k}(i,j)));
    end
    if Cstore{1,k}(1,j) < 0
        imagesc(-randomEffect)
    else
        imagesc(randomEffect)
    end
    set(gca,'YDir','normal')
    xticks(1:16)
    xticklabels(wrange)
    colormap gray
    colorbar
end

% Compared to:
for i=1:nOut
    subplot(2,nOut,i+3)
    imagesc(LMSI_real_slim(i,:))
    set(gca,'YDir','normal')
    title(labels{i})
    xticks(1:16)
    xticklabels(wrange)
    colormap gray
    colorbar
end

%%
%figure, scatter(randomEffect(1,:),LMSI_real_slim(1,:))

%%

figure,
subplot(1,nOut,1)
plot

