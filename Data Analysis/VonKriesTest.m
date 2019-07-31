% I've grabbed some code from Analyse results by cone excitation vs wavelength.m
% I want to get it to run with demo data, to show what the case would be if
% it was a simple Von Kries scaling.

%% Pre-flight

clear, clc, close all

% ciedir = fullfile('C:\Users\cege-user\Dropbox\UCL\Data\Colour standards','CIE colorimetric data');
% cvrldir = fullfile('C:\Users\cege-user\Dropbox\UCL\Data\Colour standards','CVRL cone fundamentals');

% s1 = 341;               % number of values in 1nm spectrum 390-730 nm
% w1 = 390:1:730;
% s5 = 69;                % number of values in 5nm spectrum 390-730 nm
% w5 = 390:5:730;

% Start pure --------------------

% %% Read Stockman-Sharpe 10-deg cone fundamentals (390-730 nm in 1nm intervals)
% 
% ssfile = fullfile(cvrldir,'Stockman-Sharpe cone fundamentals - lin-2deg-1nm.txt');
% format = '%d %f %f %f';
% fid = fopen(ssfile,'r');
% [Obs,count] = fscanf(fid,format,[4,inf]);  % read the whole file into array
% fclose(fid);
% 
% Lcone = Obs(2,1:s1);               % extract cone response data fields
% Mcone = Obs(3,1:s1);               % data range 380-730nm
% Scone = Obs(4,1:s1);
% LMScone = [Lcone; Mcone; Scone];         % 3x341 array
% 
% % Read scotopic CIE V'(lambda) 380-780 in 5nm intervals
% 
% vsfile = fullfile(ciedir,'CIE 1951 scotopic luminous efficiency.txt');
% format = '%d %f';
% fid = fopen(vsfile,'r');
% [V,count] = fscanf(fid,format,[2,inf]);  % read the whole file into array
% fclose(fid);
% 
% Vint = interp1(380:5:780,V(2,:),380:780,'spline');     % interpolate to 1nm
% Vprime = Vint(11:s1+10);                % extract range 390-730 nm
% 
% % End pure ---------------------

load T_cones_ss10.mat
load T_rods.mat

%% Grabbing bits to make an empty data container

wmin = 400;  wmax = 700;            % range of wavelengths (20nm intervals)
wrange = wmin:20:wmax;

N = 10;                             % number of repetitions over time
LN = 16;                            % number of lightness levels per repeat
WN = length(wrange);                % number of wavelength samples

Lval = 85:-5:10;

%TNM = floor(min(TI(4,:)/60));      % find length of shortest session (minutes)
TNM = 72;                           % DG: This is what it ends up as

LMSR = zeros(4,LN,TNM,WN,'double'); % LMSR of match (Should be 4,16,72,16)

% SO, what do these number mean?
% 4 - Long, Medium, Short, Rod
% 16 - number of lightness levels (LN)
% 72 - minutes of session
% 16 - number of wavelength samples

%% Generate simulated data

% Load spectrums of peripheral stimuli
load('C:\Users\cege-user\Dropbox\UCL\Data\LargeSphere\Hardware Data\Filter spectra\Illumination in sphere.mat','spectra')
S_spectra = [380,4,101];

for ii = 1:size(LMSR,1)         % each sensor
    for jj = 1%:LN              % (let's just do one lightness level for now)
        for kk = 1%:TNM         % (let's assume the same across time for now)
            for mm = 1:WN       % for each wavelength sample

                light = spectra(:,mm+1);
 %               sensor = SplineCmf(Obs(ii,:);
                stimulation = light * sensor;
                nsfwp = norm * stimulation; %needed stimulation for white point
                LMSR(ii,jj,kk,mm) = nsfwp;
            end
        end
    end
end


%%

% Start pure ----------------------------

%% Plot 3D surfaces of LMSR vs wavelength and time for constant lightness

lness = 8;              % lightness (fixed)
WI = 400:5:700;
wl = length(WI);
LMSI = zeros(4,wl,TNM,'double');

for t = 1:TNM
  for k = 1:4
    v = squeeze(LMSR(k,lness,t,:));    % extract LMSR values for all wavelengths
    LMSI(k,:,t) = interp1(wrange,v,WI,'spline');  % interpolate to 5nm intervals
  end
end

lp = reshape(LMSR(1,lness,:,:),TNM*WN,1);    % extract L as vector
ml = mean(lp);                              % mean L over plane
sl = std(lp);
lm = reshape(LMSR(2,lness,:,:),TNM*WN,1);    % extract M as vector
mm = mean(lm);                              % mean M over plane
sm = std(lm);
ls = reshape(LMSR(3,lness,:,:),TNM*WN,1);    % extract S as vector
ms = mean(ls);                              % mean S over plane
ss = std(ls);
fprintf(1,'\nStd and mean cone excitations for target lightness %d\n',Lval(lness));
fprintf(1,'L %f/%f=%5.3f,  M %f/%f=%5.3f,  S %f/%f=%5.3f\n',...
            sl,ml,sl/ml,sm,mm,sm/mm,ss,ms,ss/ms); 

% Plot figures

[XL,YL] = meshgrid(WI,1:TNM);

figure;  hold on;  rotate3d;  grid on;
title(sprintf('L cone excitation for wavelength vs time, lightness = %d',Lval(lness)));
ZL = squeeze(LMSI(1,:,:))';
surf(XL,YL,ZL);
xlabel('Wavelength of adapting field (nm)');
ylabel('Time (min)');
zlabel('L cone excitation');
axis([WI(1) WI(wl) 1 TNM]);
    
figure;  hold on;  rotate3d;  grid on;
title(sprintf('M cone excitation for wavelength vs time, lightness = %d',Lval(lness)));
ZL = squeeze(LMSI(2,:,:))';
surf(XL,YL,ZL);
xlabel('Wavelength of adapting field (nm)');
ylabel('Time (min)');
zlabel('M cone excitation');
axis([WI(1) WI(wl) 1 TNM]);

figure;  hold on;  rotate3d;  grid on;
title(sprintf('S cone excitation for wavelength vs time, lightness = %d',Lval(lness)));
ZL = squeeze(LMSI(3,:,:))';
surf(XL,YL,ZL);
xlabel('Wavelength of adapting field (nm)');
ylabel('Time (min)');
zlabel('S cone excitation');
axis([WI(1) WI(wl) 1 TNM]);

figure;  hold on;  rotate3d;  grid on;
title(sprintf('Rod excitation for wavelength vs time, lightness = %d',Lval(lness)));
ZL = squeeze(LMSI(4,:,:))';
surf(XL,YL,ZL);
xlabel('Wavelength of adapting field (nm)');
ylabel('Time (min)');
zlabel('Rod excitation');
axis([WI(1) WI(wl) 1 TNM]);
