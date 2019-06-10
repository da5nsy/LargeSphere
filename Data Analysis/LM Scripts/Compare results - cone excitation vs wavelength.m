%--------------------------------------------------------------------------
%  VISION SPHERE EXPERIMENT - WHITE BALANCE - TIME SERIES
%
%  Compare two sets of experimental results for
%  cone excitations over all wavelengths
%
%--------------------------------------------------------------------------

load('Lindsay LMS');                % load experimental results #1
LABM1 = LABM;
RGBM1 = RGBM;
TM1 = TM;
RGBL1 = RGBL;
RGBlut1 = RGBlut;
TI1 = TI;
LMSR1 = LMSR;
TNM1 = floor(min(TI1(4,:)/60));     % length of shortest session (minutes)

load('Tania LMS');                % load experimental results #2
LABM2 = LABM;
RGBM2 = RGBM;
TM2 = TM;
RGBL2 = RGBL;
RGBlut2 = RGBlut;
TI2 = TI;
LMSR2 = LMSR;
TNM2 = floor(min(TI2(4,:)/60));     % length of shortest session (minutes)

wmin = 400;  wmax = 700;            % range of wavelengths (20nm intervals)
wrange = wmin:20:wmax;
N = 10;                             % number of repetitions over time
LN = 16;                            % number of lightness levels per repeat
WN = length(wrange);                % number of wavelength samples

%% Plot session lengths

min1 = min(TI1(4,:))/60;  max1 = max(TI1(4,:))/60;
min2 = min(TI2(4,:))/60;  max2 = max(TI2(4,:))/60;
fprintf('Range of session lengths #1 %d-%d  #2 %d-%d min\n',...
    floor(min1),floor(max1),floor(min2),floor(max2));
ms1 = sum(TI1(4,:)/(16*60));  ms2 = sum(TI2(4,:)/(16*60));
fprintf('Mean session length #1 %d  #2 %d min\n',floor(ms1),floor(ms2));

figure;  hold on;
axis([400 700 0 120]);
plot(400:20:700,TI2(4,:)/60,'o-m');
plot(400:20:700,TI1(4,:)/60,'o-b');
legend('TR','LM');
plot([400 700],[ms1 ms1],':b');
plot([400 700],[ms2 ms2],':m');
xlabel('Adapting wavelength (nm)');
ylabel('Session length (mins)');

%% Plot LMS over time for one wavelength

lam = 3;
lambda = wrange(lam);

figure;  hold on;
title(sprintf('LMS cone excitation vs time, adapting wavelength = %d nm',lambda));

% #1

tmax1 = TI1(4,lam);                       % session length (sec) for this wavelength
xlim([1 TNM1]);                           % X axis in minutes
rvi = squeeze(LMSR1(1,:,1:TNM1,lam));     % extract all L values for wavelength
rvm = mean(rvi(1:LN,:));                  % mean over all lightnesses
plot(1:TNM1,rvm,'--r','LineWidth',2);      % plot as thick line
gvi = squeeze(LMSR1(2,:,1:TNM1,lam));     % extract all M values for wavelength
gvm = mean(gvi(1:LN,:));                  % mean over all lightnesses
plot(1:TNM1,gvm,'--g','LineWidth',2);      % plot as thick line
bvi = squeeze(LMSR1(3,:,1:TNM1,lam));     % extract all S values for wavelength
bvm = mean(bvi(1:LN,:));                  % mean over all lightnesses
plot(1:TNM1,bvm,'--b','LineWidth',2);      % plot as thick line

% #2

tmax2 = TI2(4,lam);                       % session length (sec) for this wavelength
xlim([1 TNM2]);                           % X axis in minutes
rvi = squeeze(LMSR2(1,:,1:TNM2,lam));     % extract all L values for wavelength
rvm = mean(rvi(1:LN,:));                  % mean over all lightnesses
plot(1:TNM2,rvm,'-r','LineWidth',2);      % plot as thick line
gvi = squeeze(LMSR2(2,:,1:TNM2,lam));     % extract all M values for wavelength
gvm = mean(gvi(1:LN,:));                  % mean over all lightnesses
plot(1:TNM2,gvm,'-g','LineWidth',2);      % plot as thick line
bvi = squeeze(LMSR2(3,:,1:TNM2,lam));     % extract all S values for wavelength
bvm = mean(bvi(1:LN,:));                  % mean over all lightnesses
plot(1:TNM2,bvm,'-b','LineWidth',2);      % plot as thick line

xlabel('Time (min)');
ylabel('Value relative to mean over all iterations');
