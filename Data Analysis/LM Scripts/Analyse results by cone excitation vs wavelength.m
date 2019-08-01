%-------------------------------------------------------------------------
%  VISION SPHERE EXPERIMENT - WHITE BALANCE - TIME SERIES
%
%  Analyse experimental results in terms of L,M,S cone excitations
%
%-------------------------------------------------------------------------

clear, clc, close all

rootdir = 'C:\Users\cege-user\Dropbox\UCL\Data\LargeSphere';
dir = 'C:\Users\cege-user\Dropbox\UCL\Data\LargeSphere\Experimental Data\2013 Apr TR';
%dir = 'C:\Users\cege-user\Dropbox\UCL\Data\LargeSphere\Experimental Data\2012 Sep LM\Lindsay time series - Sep 2012';

s1 = 341;               % number of values in 1nm spectrum 390-730 nm
w1 = 390:1:730;

% Load observer
load T_cones_ss10.mat
load T_rods.mat
T_obs = [T_cones_ss10; SplineCmf(S_rods,T_rods,S_cones_ss10)];
S_obs = S_cones_ss10;

LMScone = T_cones_ss10(:,1:341);
Vprime = SplineCmf(S_rods,T_rods,w1'); %for back compatibility with LM's code
V = [SToWls(S_rods)';T_rods]; %for back compatibility with LM's code

clear T_cones_ss10 S_cones_ss10 T_rods S_rods %cleanup
figure, plot(SToWls(S_obs),T_obs')


%% Read experimental data for neutral matching

wmin = 400;  wmax = 700;            % range of wavelengths (20nm intervals)
wrange = wmin:20:wmax;
N = 10;                             % number of repetitions over time
LN = 16;                            % number of lightness levels per repeat
WN = length(wrange);                % number of wavelength samples

LABM = zeros(3,LN,N,WN,'double');
RGBM = zeros(3,LN,N,WN,'double');
TM = zeros(3,LN,N,WN,'double');

for n = 1:WN
  lam = wrange(n);                  % filter wavelength
  fname = fullfile(dir,sprintf('%dnm - time',lam));
  load(fname);                      % load experimental results for wavelength
  LABM(:,:,:,n) = LABmatch(:,1:LN,:);
  RGBM(:,:,:,n) = RGBmatch(:,1:LN,:);   % copy all data into a single array
  TM(:,:,:,n) = Tmatch(:,1:LN,:);
end

Lval = squeeze(LABM(1,:,1,1));      % L values

%% Read display data for primaries and luminance calibation

fname = 'C:\Users\cege-user\Dropbox\UCL\Data\LargeSphere\Hardware Data\Calibration Data\Large LCD display measurement.mat';
load(fname);                        % read data from file

imax = 65535;                       % maximum integer in 16-bit range
dac = 0:imax;                       % range of input values
ndac = dac/imax;                    % normalised to range 0-1
RGBlut = zeros(3,imax+1,'double');  % lookup tables for tone curve (quasi 16-bit)
RGBL = zeros(3,LN,N,WN,'double');   % converted normalised luminance values

for k = 1:3
  Yval = squeeze(XYZ(2,:,k));       % get vector of Y values
  RGBlut(k,:) = interp1(256*sval,Yval,dac,'spline');  % interpolate Y values
end

Rlut = squeeze(RGBlut(1,:));        % extract 1-D versions
Glut = squeeze(RGBlut(2,:));
Blut = squeeze(RGBlut(3,:));

rw = RGBW_XYZ(2,1);                 % luminance of R,G,B display primaries (max signal)
gw = RGBW_XYZ(2,2);
bw = RGBW_XYZ(2,3);

%% Convert RGB match values to display luminance

for n = 1:WN
 for t = 1:N
  for i = 1:LN
    r = RGBM(1,i,t,n);
    ri = 1+uint16(65534*r);
    RGBL(1,i,t,n) = Rlut(ri)/rw;   % convert to fraction of white luminance
    g = RGBM(2,i,t,n);
    gi = 1+uint16(65534*g);
    RGBL(2,i,t,n) = Glut(gi)/gw;
    b = RGBM(3,i,t,n);
    bi = 1+uint16(65534*b);
    RGBL(3,i,t,n) = Blut(bi)/bw;
  end
 end
end

%% Interpolate display primary spectra to 1nm intervals

DS = zeros(s1,3,'double');          % RGB display spectra
w4 = lambda;                        % measurement wavelengths (380-780, 4nm intervals)

for k = 1:3
  p = squeeze(RGBW_spectrum(:,k));
  DS(:,k) = interp1(w4,p,w1','spline');
end

figure;  hold on;
title('Display RGB spectra');
plot(w1,DS(:,1),'-r');
plot(w1,DS(:,2),'-g');
plot(w1,DS(:,3),'-b');
xlabel('Wavelength (nm)');
ylabel('Power');

%% Analyse relative RGB vs time

TI = zeros(4,WN,'double');          % start time (h,m,s) and elapsed time (sec)

% Calculate session lengths for each wavelength

for n = 1:WN
  t1h = squeeze(TM(1,1,1,n));               % time of first observation
  t1m = squeeze(TM(2,1,1,n));
  t1s = squeeze(TM(3,1,1,n));
  t2h = squeeze(TM(1,LN,N,n));              % time of last observation
  t2m = squeeze(TM(2,LN,N,n));
  t2s = squeeze(TM(3,LN,N,n));
  if n == 14
    t2h = squeeze(TM(1,LN-1,N,n));          % special case for 660 nm
    t2m = squeeze(TM(2,LN-1,N,n));          % (missing lowest lightness level)
    t2s = squeeze(TM(3,LN-1,N,n));
  end
  TI(1,n) = t1h;
  TI(2,n) = t1m;
  TI(3,n) = t1s;
  TI(4,n) = round(3600*(t2h-t1h)+60*(t2m-t1m)+(t2s-t1s));   % length of session from first to last (sec)
end

fprintf('Range of session lengths %d-%d min\n',floor(min(TI(4,:)/60)),floor(max(TI(4,:)/60)));
ms = sum(TI(4,:)/(16*60));
fprintf('Mean session length %d min\n',floor(ms));

figure;  hold on;
axis([400 700 0 120]);
plot([400 700],[ms ms],':k');
plot(400:20:700,TI(4,:)/60,'o-b');
xlabel('Adapting wavelength (nm)');
ylabel('Session length (mins)');

%% Interpolate RGB to relative value every second

TNS = max(TI(4,:));                        % maximum session length (sec)
RGBval = zeros(3,LN,TNS,WN,'double');      % actual value per second (range 0-1)
RGBrel = zeros(3,LN,TNS,WN,'double');      % value relative to mean over iterations

for n = 1:WN
 t1h = TI(1,n);
 t1m = TI(2,n);
 t1s = TI(3,n);
 tmax = TI(4,n);
 for k = 1:3
  for i = 1:LN
   v = squeeze(RGBM(k,i,:,n));      % match values for one lightness, all iterations
   th = squeeze(TM(1,i,:,n));       % hours
   tm = squeeze(TM(2,i,:,n));       % minutes
   ts = squeeze(TM(3,i,:,n));       % seconds
   te = 3600*(th-t1h) + 60*(tm-t1m) + (ts-t1s);  % elapsed time for each observation (sec)
   vi = interp1(te,v,1:tmax,'linear','extrap');  % interpolate to seconds  
   RGBval(k,i,1:tmax,n) = vi;       % save interpolated lightness values
   vv = squeeze(RGBM(k,i,:,:));     % match values for one lightness, all iterations, all wavelengths
   T = vv<0;  vv(T) = 0;            % clip negative values to zero
   mv = geomean(geomean(vv));       % mean  over all iterations and wavelengths
   RGBrel(k,i,1:tmax,n) = vi/mv;    % normalise to average over all iterations
  end
 end
end

%% Convert display RGB to true tristimulus XYZ at 1 second intervals

Xw65 = 95.04; Yw65 = 100; Zw65 = 108.86;  % white tristimulus for D65
XwD = 99.04; YwD = 100; ZwD = 151.30;     % normalised white tristimulus for display

rw = RGBW_XYZ(2,1);     % luminance of R,G,B display primaries (max signal)
gw = RGBW_XYZ(2,2);
bw = RGBW_XYZ(2,3);
M = RGBW_XYZ(1:3,1:3);  % matrix of XYZ vs RGB

TNM = floor(min(TI(4,:)/60));           % find length of shortest session (minutes)
XYZD = zeros(3,LN,TNM,WN,'double');     % display XYZ of match
XYZD65 = zeros(3,LN,TNM,WN,'double');   % display XYZ of match (D65 white ref)
LABD = zeros(3,LN,TNM,WN,'double');     % display LAB of match
RGBD = zeros(3,LN,TNM,WN,'uint8');      % display RGB of match (sRGB, D65)

for n = 1:WN
 for t = 1:TNM
  for i = 1:LN
    ts = 60*t;                  % time in seconds
    rgb = RGBval(:,i,ts,n);     % get RGB value (display signal, range 0-1)
    rgbi = 1+uint16(65534*rgb); % convert to 16-bit LUT index
    rl = Rlut(rgbi(1))/rw;      % convert to fraction of white luminance
    gl = Glut(rgbi(2))/gw;
    bl = Blut(rgbi(3))/bw;
    xyz = M*[rl gl bl]';        % tristimulus values
    XYZD(:,i,t,n) = xyz;
    X = xyz(1);
    Y = xyz(2);
    Z = xyz(3);
    [L A B] = XYZtoLAB(X,Y,Z,XwD,YwD,ZwD);          % convert to LAB with display white
    LABD(:,i,t,n) = [L A B];
    [X65 Y65 Z65] = LABtoXYZ(L,A,B,Xw65,Yw65,Zw65); % convert for D65 white ref
    XYZD65(:,i,t,n) = [X65 Y65 Z65];
    [r g b gamutflag] = XYZtosRGBe(X65/100,Y65/100,Z65/100); % convert to sRGB
    RGBD(:,i,t,n) = uint8(255*[r g b]);
  end
 end
end

%% Calculate LMSR response to RGB display primaries

RGB_LMSR = zeros(4,3,'double');  % RGB columns, LMSR rows

for i = 1:3
  for k = 1:3
    RGB_LMSR(k,i) = sum(DS(:,i)'.*LMScone(k,:));  % sum over range 390-730nm
  end
end

for i = 1:3
  RGB_LMSR(4,i) = sum(DS(:,i)'.*Vprime);  % rod response
end

RGB_LMSR = RGB_LMSR/RGB_LMSR(2,2);     % normalise to green on M cone

disp(RGB_LMSR);

%% Convert display RGB to LMSR excitations at 1 minute intervals

rgbw = RGBW_XYZ(2,:);     % luminance of R,G,B display primaries (max signal)

TNM = floor(min(TI(4,:)/60));       % find length of shortest session (minutes)
LMSR = zeros(4,LN,TNM,WN,'double');  % LMSR of match

for n = 1:WN
 for t = 1:TNM
  for i = 1:LN
    ts = 60*t;                  % time in seconds
    rgb = RGBval(:,i,ts,n);     % get RGB match value (display signal, range 0-1)
    rgbi = 1+uint16(65534*rgb); % convert to 16-bit LUT index
    rl = Rlut(rgbi(1))/rw;      % convert to fraction of white luminance
    gl = Glut(rgbi(2))/gw;
    bl = Blut(rgbi(3))/bw;
    rgbl = [rl gl bl];          % RGB luminances
    for k = 1:4
      lmsrvec = RGB_LMSR(k,:)';   % row of matrix for one of L,M,S,R
      LMSR(k,i,t,n) = rgbl*lmsrvec;   % sum contribution from each primary
    end
  end
 end
end

%save('Lindsay LMS','LMSR','LABM','RGBM','TM','TI','RGBlut','RGBL');

%% Plot cone excitations in LMS cube

p = 1;                              % position of axis markers
m = 10;                             % size of axis markers

lmsmax = max(max(max(max(LMSR(1:3,:,:,:)))));   % max cone value for scaling

figure;  hold on;  grid on;  rotate3d on;  axis equal;
title('LMS cone excitation values for all matches');

for w = 1:WN
  for n = 1:N
    for i = 1:LN
       lms = LMSR(:,i,n,w)/lmsmax;        % normalised cone values (range 0-1)
       rgb = double(RGBD(:,i,n,w))/255;  % sRGB display values (range 0-1)
       T = rgb<0;  rgb(T) = 0;           % clip to zero
       plot3(lms(1),lms(2),lms(3),'o','MarkerSize',5,'MarkerEdgeColor',rgb,...
      'MarkerFaceColor',rgb);
    end
  end
end

axis([0 p 0 p 0 p]);
plot3([0 p],[0 p],[0 p],'--k');  % draw diagonal of cube
plot3([0 p],[0 0],[0 0],'-k');   % draw primary edges of cube
plot3([0 0],[0 p],[0 0],'-k');
plot3([0 0],[0 0],[0 p],'-k');
plot3([0 p],[0 0],[p p],'-k');   % draw outer edges of cube
plot3([p p],[0 p],[0 0],'-k');
plot3([0 0],[p p],[0 p],'-k');
plot3([0 p],[p p],[0 0],'-k');   % draw inner edges of cube
plot3([0 0],[0 p],[p p],'-k');
plot3([p p],[0 0],[0 p],'-k');
plot3([0 p],[p p],[p p],'-k');   % draw outer edges of cube
plot3([p p],[0 p],[p p],'-k');
plot3([p p],[p p],[0 p],'-k');
plot3(0,0,0,'ok','MarkerSize',m,'MarkerFaceColor','k');
plot3(p,p,p,'ok','MarkerSize',m,'MarkerFaceColor','w');
xlabel('L cone');
ylabel('M cone');
zlabel('S cone');

%% Plot 3D surfaces of LMSR vs wavelength and time for constant lightness

lness = 6:11;              % lightness (can be single or multiple)
LMSI=squeeze(mean(LMSR(:,lness,:,:),2));
labels = {'L','M','S'};

figure('units','normalized','outerposition',[0 0 1 1])
for i=1:3
    subplot(1,3,i)
    imagesc(squeeze(LMSI(i,:,:)))
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

%% Plot figures

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

% %% Plot 3D surfaces vs wavelength and lightness at given time
% 
% tm = 30;                              % time (minutes)
% LNI = 20:75;
% li = length(LNI);
% LMSI = zeros(3,wl,li,'double');
% 
% for k = 1:3
%   V = squeeze(LMSR(k,3:LN-2,tm,:));    % extract LMS values for all lightness levels and wavelengths
%   for i = 1:LN-4
%     vm = mean(V(i,:));                % mean value over all wavelengths for each lightness
%     V(i,:) = V(i,:)/vm;
%   end
%   LMSI(k,:,:) = interp2(Lval(3:LN-2)',wrange,V',LNI',WI,'spline');  % interpolate to 5nm intervals
% end
% 
% [XL,YL] = meshgrid(WI,LNI);
% 
% figure;  hold on;  rotate3d;  grid on;
% title(sprintf('L cone excitation for wavelength vs lightness, time = %d min',tm));
% ZL = squeeze(LMSI(1,:,:))';
% surf(XL,YL,ZL);
% axis([400 700 20 75 0.8 1.2]);
% xlabel('Wavelength of adapting field (nm)');
% ylabel('Lightness of target');
% zlabel('L cone excitation relative to mean across all wavelengths');
% 
% figure;  hold on;  rotate3d;  grid on;
% title(sprintf('M cone excitation for wavelength vs lightness, time = %d min',tm));
% ZL = squeeze(LMSI(2,:,:))';
% surf(XL,YL,ZL);
% axis([400 700 20 75 0.8 1.2]);
% xlabel('Wavelength of adapting field (nm)');
% ylabel('Lightness of target');
% zlabel('M cone excitation relative to mean across all wavelengths');
% 
% figure;  hold on;  rotate3d;  grid on;
% title(sprintf('S cone excitation for wavelength vs lightness, time = %d min',tm));
% ZL = squeeze(LMSI(3,:,:))';
% surf(XL,YL,ZL);
% axis([400 700 20 75 0.5 1.5]);
% xlabel('Wavelength of adapting field (nm)');
% ylabel('Lightness of target');
% zlabel('S cone excitation relative to mean across all wavelengths');
% 
% %% Plot surface of L-M chromatic opponent channel for wavelength vs time
% 
% lness = 8;
% YCC = zeros(3,wl,TNM,'double');   % Y=R+G, Cr=R-G, Cb=Y-B
% 
% for t = 1:TNM
%   cr = squeeze(LMSR(1,lness,t,:));    % extract L values for all wavelengths
%   cri = interp1(wrange,cr,WI,'spline');
%   cg = squeeze(LMSR(2,lness,t,:));    % extract M values for all wavelengths
%   cgi = interp1(wrange,cg,WI,'spline');
%   cb = squeeze(LMSR(3,lness,t,:));    % extract S values for all wavelengths
%   cbi = interp1(wrange,cb,WI,'spline');
%   YCC(1,:,t) = cri+cgi;
%   YCC(2,:,t) = cri-cgi;
%   YCC(3,:,t) = cri+cgi-cbi;
% end
% 
% [XL,YL] = meshgrid(WI,1:TNM);
% 
% figure;  hold on;  rotate3d;  grid on;
% ylim([1 TNM]);
% title(sprintf('Cone L+M luminance for wavelength vs time, lightness = %d',Lval(lness)));
% ZL = squeeze(YCC(1,:,:));
% surf(XL,YL,ZL');
% xlabel('Wavelength of adapting field (nm)');
% ylabel('Time (min)');
% zlabel('Cone L+M');
% 
% figure;  hold on;  rotate3d;  grid on;
% title(sprintf('Cone L-M opponent channel for wavelength vs time, lightness = %d',Lval(lness)));
% ylim([1 TNM]);
% ZL = squeeze(YCC(2,:,:));
% surf(XL,YL,ZL');
% xlabel('Wavelength of adapting field (nm)');
% ylabel('Time (min)');
% zlabel('Cone L-M');
% 
% figure;  hold on;  rotate3d;  grid on;
% title(sprintf('Cone Y-S opponent channel for wavelength vs time, lightness = %d',Lval(lness)));
% ylim([1 TNM]);
% ZL = squeeze(YCC(3,:,:));
% surf(XL,YL,ZL');
% xlabel('Wavelength of adapting field (nm)');
% ylabel('Time (min)');
% zlabel('Cone L+M-S');
% 
% %% Plot LMS over time for one wavelength
% 
% lam = 2;
% lambda = wrange(lam);
% 
% figure;  hold on;
% title(sprintf('LMS cone excitation vs time, adapting wavelength = %d nm',lambda));
% tmax = TI(4,lam);                         % session length (sec) for this wavelength
% xlim([1 TNM]);                            % X axis in minutes
% rvi = squeeze(LMSR(1,:,1:TNM,lam));       % extract all L values for wavelength
% rvm = mean(rvi(1:LN,:));                  % mean over all lightnesses
% plot(1:TNM,rvm,'-r','LineWidth',2);       % plot as thick line
% gvi = squeeze(LMSR(2,:,1:TNM,lam));       % extract all M values for wavelength
% gvm = mean(gvi(1:LN,:));                  % mean over all lightnesses
% plot(1:TNM,gvm,'-g','LineWidth',2);       % plot as thick line
% bvi = squeeze(LMSR(3,:,1:TNM,lam));       % extract all S values for wavelength
% bvm = mean(bvi(1:LN,:));                  % mean over all lightnesses
% plot(1:TNM,bvm,'-b','LineWidth',2);       % plot as thick line
% %kvi = squeeze(LMSR(4,:,1:TNM,lam));      % extract all R values for wavelength
% %kvm = mean(kvi(1:LN,:));                 % mean over all lightnesses
% %plot(1:TNM,kvm,'-k','LineWidth',2);      % plot as thick line
% tm = ceil(max(TI(4,:))/60);               % length of longest session (minutes)
% xlabel('Time (min)');
% ylabel('Value relative to mean over all iterations');
