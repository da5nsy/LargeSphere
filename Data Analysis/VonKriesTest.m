% I've grabbed some code from Analyse results by cone excitation vs wavelength.m
% I want to get it to run with demo data, to show what the case would be if
% it was a simple Von Kries scaling.

%%

clear, clc, close all

plt = 0;

% Load observer
load T_cones_ss10.mat
load T_rods.mat
T_obs = [T_cones_ss10; SplineCmf(S_rods,T_rods,S_cones_ss10)];
S_obs = S_cones_ss10;
clear T_cones_ss10 S_cones_ss10 T_rods S_rods %cleanup
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

LMSR_sim = zeros(4,LN,TNM,WN,'double'); % Empty data container (Should be 4,16,72,16)

% SO, what do these number mean?

% 4 - Long, Medium, Short, Rod
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
labels = {'L','M','S'};

if plt
    figure('units','normalized','outerposition',[0 0 1 1])
    for i=1:3
        subplot(1,3,i)
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

load 'Tania LMSI for testing VonKriesTest'
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

%rng(1)

storeCW_L = zeros(4,1);
storeCW_M = zeros(4,1);
storeCW_S = zeros(4,1);
storeCW_R = zeros(4,1);

for j=1:10000
    ConeWeights = rand(size(LMSI_sim,1),1)*20-10;
    randomEffect = zeros([72,16]);
    for i=1:3
        randomEffect = randomEffect + squeeze((LMSI_sim(i,:,:)*ConeWeights(i,1)));
    end
    
%     figure
%     imagesc(randomEffect)
%     set(gca,'YDir','normal')
%     
%     xticks(1:16)
%     xticklabels(wrange)
%     colormap gray
    
    randomEffect_slim = mean(randomEffect,1);
    % Compare random effect with our actual data
    LMSI_real_slim = squeeze(mean(LMSI_real,2));
    for k=1:3
        if min(corrcoef(randomEffect_slim,LMSI_real_slim(k,:)')) >0.77
            corrcoef(randomEffect_slim,LMSI_real_slim(k,:)');
            if k==1
            storeCW_L = cat(3,storeCW_L,ConeWeights);
            elseif k==2                
            storeCW_M = cat(3,storeCW_M,ConeWeights);
            elseif k==3
            storeCW_S = cat(3,storeCW_S,ConeWeights);
%             elseif k==4
%             storeCW_R = cat(3,storeCW_R,ConeWeights);
            end
        end
        
    end
end

%% Visualise successes
for j=2%:10%size(storeCW,3)
    randomEffect = zeros([72,16]);
    for i=1:3
        randomEffect = randomEffect + squeeze((LMSI_sim(i,:,:)*storeCW_L(i,1,j)));
    end
    figure
    imagesc(randomEffect)
    set(gca,'YDir','normal')
    
    xticks(1:16)
    xticklabels(wrange)
    colormap gray
    colorbar
end

for j=2%:10%size(storeCW,3)
    randomEffect = zeros([72,16]);
    for i=1:3
        randomEffect = randomEffect + squeeze((LMSI_sim(i,:,:)*storeCW_M(i,1,j)));
    end
    figure
    imagesc(randomEffect)
    set(gca,'YDir','normal')
    
    xticks(1:16)
    xticklabels(wrange)
    colormap gray
    colorbar
end

for j=2%:10%size(storeCW,3)
    randomEffect = zeros([72,16]);
    for i=1:3
        randomEffect = randomEffect + squeeze((LMSI_sim(i,:,:)*storeCW_S(i,1,j)));
    end
    figure
    imagesc(randomEffect)
    set(gca,'YDir','normal')
    
    xticks(1:16)
    xticklabels(wrange)
    colormap gray
    colorbar
end

%% Compared to:
figure,
imagesc(LMSI_real_slim(1,:))
set(gca,'YDir','normal')
title('REAL(ish) L')
xticks(1:16)
xticklabels(wrange)
colormap gray

figure,
imagesc(LMSI_real_slim(2,:))
set(gca,'YDir','normal')
title('REAL(ish) M')
xticks(1:16)
xticklabels(wrange)
colormap gray

figure,
imagesc(LMSI_real_slim(3,:))
set(gca,'YDir','normal')
title('REAL(ish) S')
xticks(1:16)
xticklabels(wrange)
colormap gray


