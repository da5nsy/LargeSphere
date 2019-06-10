%-------------------------------------------------------------------------
%  VISION SPHERE EXPERIMENT - WHITE BALANCE - TIME SERIES
%
%  Analyse experimental results for one wavelength and plot figures
%
%-------------------------------------------------------------------------

clear, clc, close all

name = 'Lindsay';
filter_lambda = 480;                  % filter wavelength

rootdir = fullfile('C:','Users','cege-user','Dropbox','UCL','Data','LargeSphere');

if strcmp(name,'Lindsay');
  fname = sprintf('%dnm - time.mat',filter_lambda);
  ifile = fullfile(rootdir,'Lindsay time series - Sep 2012',fname);
elseif strcmp(name,'Tania');
  fname = sprintf('%dnm - time.mat',filter_lambda);
  ifile = fullfile(rootdir,'Tania time series - Apr 2013',fname);
else
  fname = sprintf('%dnm - time v2_4.mat',filter_lambda);
  ifile = fullfile(rootdir,'Danny time series - Dec 2016','Results',fname);
end
load(ifile);                 % load experimental results for wavelength

N = 10;                             % number of repetitions over time
LN = 16;                            % number of lightness levels per repeat

Lval = squeeze(LABmatch(1,:,1));    % L values
dfile = fullfile(rootdir,'Large LCD display measurement.mat');
load(dfile);                        % load display data

%% Make mosaic of matching colours - lightness vs iteration

b = 40;                         % pixels in box side
s = 4;                          % spacing between boxes
w = s+N*(b+s);                  % width of array (iteration axis)
h = s+LN*(b+s);                 % height of array (lightness axis)
Im = zeros(h,w,3,'uint8');      % image array

for n = 1:N
  xp = s+(n-1)*(b+s);                    % x pixel address (iteration axis)
  for i = 1:LN
    rs = RGBmatch(1,i,n);                % get R value (display signal, 8-bit)
    gs = RGBmatch(2,i,n);
    bs = RGBmatch(3,i,n);
    yp = s+(i-1)*(b+s);                  % y pixel address (lightness axis)
    Im(yp:yp+b-1,xp:xp+b-1,1) = uint8(255*rs);  % fill one square in array
    Im(yp:yp+b-1,xp:xp+b-1,2) = uint8(255*gs);
    Im(yp:yp+b-1,xp:xp+b-1,3) = uint8(255*bs);
  end
end

iname = fullfile(rootdir,sprintf('%s matches, %dnm_4.tif',name,filter_lambda));
imwrite(Im,iname,'tif');           % write the image

%% Plot results as 3D surfaces

LNR = 2:LN-2;                       % range of lightness 80-20
figure;  hold on;  rotate3d;  grid on;
[XL,YL] = meshgrid(Lval(LNR),1:10); % lightness on X axis; iteration on Y axis
ZL = squeeze(RGBmatch(1,LNR,:));    % red display values
T = ZL<0;  ZL(T) = 0;               % clip negative values to zero
for n = LNR-1
  rm = mean(ZL(n,:));
  ZL(n,:) = ZL(n,:)/rm;             % normalise over mean of all iterations
end
surf(XL,YL,ZL');
title(sprintf('%s: %s, Red display signal',name,fname));
xlabel('Lightness');
ylabel('Iteration');
zlabel('Red display signal relative to mean');
axis([Lval(LN-2) Lval(2) 1 N min(min(ZL)) max(max(ZL))]);

figure;  hold on;  rotate3d;  grid on;
ZL = squeeze(RGBmatch(2,LNR,:));         % green display values
T = ZL<0;  ZL(T) = 0;                  % clip negative values to zero
for n = LNR-1
  gm = mean(ZL(n,:));
  ZL(n,:) = ZL(n,:)/gm;                % normalise over mean of all iterations
end
surf(XL,YL,ZL');
title(sprintf('%s: %s, Green display signal',name,fname));
xlabel('Lightness');
ylabel('Iteration');
zlabel('Green display signal relative to mean');
axis([Lval(LN-2) Lval(2) 1 N min(min(ZL)) max(max(ZL))]);

figure;  hold on;  rotate3d;  grid on;
ZL = squeeze(RGBmatch(3,LNR,:));          % blue display values
T = ZL<0;  ZL(T) = 0;                 % clip negative values to zero
for n = LNR-1
  bm = mean(ZL(n,:));
  ZL(n,:) = ZL(n,:)/bm;                % normalise over mean of all iterations
end
surf(XL,YL,ZL');
title(sprintf('%s: %s, Blue display signal',name,fname));
xlabel('Lightness');
ylabel('Iteration');
zlabel('Blue display signal relative to mean');
axis([Lval(LN-2) Lval(2) 1 N min(min(ZL)) max(max(ZL))]);

%%  Analyse LAB differences vs iteration

figure;  hold on;
title(sprintf('A slider values vs iteration, wavelength = %d nm',filter_lambda));
for n = 1:LN
  av = squeeze(LABmatch(2,n,:));
  plot(1:N,av,'-r');
end
xlabel('Iteration');
ylabel('A');

figure;  hold on;
title(sprintf('B slider values vs iteration, wavelength = %d nm',filter_lambda));
for n = 1:LN-1
  bv = squeeze(LABmatch(3,n,:));
  plot(1:N,bv,'-b');
end
xlabel('Iteration');
ylabel('B');

%% R/B ratio vs iteration

figure;  hold on;
title(sprintf('Ratio R/B display values vs iteration, wavelength = %d nm',filter_lambda));
for n = 1:LN-2
  rv = squeeze(RGBmatch(1,n,:));
  bv = squeeze(RGBmatch(3,n,:));
  plot(1:N,rv./bv,'-k');            % plot for individual iterations
end
plot([1 N],[1 1],':k');             % dotted line at one
for n = 1:N
  rv = squeeze(RGBmatch(1,1:LN-2,n))';
  bv = squeeze(RGBmatch(3,1:LN-2,n))';
  rbm(n) = mean(rv./bv);
end
plot(1:N,rbm,'-m','LineWidth',3);   % plot mean of all iterations
xlabel('Iteration');
ylabel('R/B');

%% Plot white balance signals vs lightness

figure;  hold on;
title(sprintf('LAB slider settings for white balance, %d nm',filter_lambda));
axis([10 85 -50 50]);
for n = 1:N
  plot(LABmatch(1,:,1),LABmatch(2,:,n),'-r');
  plot(LABmatch(1,:,1),LABmatch(3,:,n),'-b');
end
xlabel('Lightness');
ylabel('Slider setting');
legend('A','B','Location','SouthWest');
plot([10 85],[0 0],':k');

% RGB

figure;  hold on;
title(sprintf('RGB display signals for white balance, %d nm',filter_lambda));
axis([10 85 0 1]);
for n = 1:N
  plot(LABmatch(1,:,1),RGBmatch(1,:,n),'-r');
  plot(LABmatch(1,:,1),RGBmatch(2,:,n),'-g');
  plot(LABmatch(1,:,1),RGBmatch(3,:,n),'-b');
end
xlabel('Lightness');
ylabel('Normalised display values');
legend('R','G','B','Location','NorthWest');

%% Analyse RGB signals vs iteration

figure;  hold on;
title(sprintf('Red display values vs iteration, wavelength = %d nm',filter_lambda));
for n = 1:LN
  rv = squeeze(RGBmatch(1,n,:));
  plot(1:N,rv,'-r');
end
xlabel('Iteration');
ylabel('R');

figure;  hold on;
title(sprintf('Green display values vs iteration, wavelength = %d nm',filter_lambda));
for n = 1:LN
  gv = squeeze(RGBmatch(2,n,:));
  plot(1:N,gv,'-g');
end
xlabel('Iteration');
ylabel('G');

figure;  hold on;
title(sprintf('Blue display values vs iteration, wavelength = %d nm',filter_lambda));
for n = 1:LN
  bv = squeeze(RGBmatch(3,n,:));
  plot(1:N,bv,'-b');
end
xlabel('Iteration');
ylabel('B');

%% Interpolate relative RGB values vs time

t1h = squeeze(Tmatch(1,1,1));               % time of first observation
t1m = squeeze(Tmatch(2,1,1));
t1s = squeeze(Tmatch(3,1,1));
t2h = squeeze(Tmatch(1,LN,N));              % time of last observation
t2m = squeeze(Tmatch(2,LN,N));
t2s = squeeze(Tmatch(3,LN,N));
tmax = round(3600*(t2h-t1h)+60*(t2m-t1m)+(t2s-t1s));   % length of session from first to last (sec)
trange = 1:tmax;
fprintf('Session length = %d minutes\n',round(tmax/60));
fprintf('Mean time per observation = %d seconds\n',round(tmax/(LN*N)));

% Red

figure;  hold on;
title(sprintf('Red display values vs time, wavelength = %d nm',filter_lambda));
RS = zeros(tmax,LN,'double');
rvi = zeros(tmax,LN,'double');

for n = 1:LN
  rv = squeeze(RGBmatch(1,n,:));       % R values for one lightness, all iterations
  rm = mean(rv);                       % mean over iterations
  th = squeeze(Tmatch(1,n,:));
  tm = squeeze(Tmatch(2,n,:));
  ts = squeeze(Tmatch(3,n,:));
  te = 3600*(th-t1h) + 60*(tm-t1m) + (ts-t1s);  % elapsed time for each observation (sec)
  RS(:,n) = interp1(te,rv,1:tmax,'linear','extrap');  % interpolate to seconds  
  rvi(:,n) = RS(:,n)'/rm;                   % normalise to average over all iterations
  plot(te/60,rv/rm,'-k');
end
rvm = mean(rvi(:,2:LN-1),2);                % mean over all lightnesses except first and last
plot(trange/60,rvm,'-r','LineWidth',3);     % plot as thick red line
plot([1 ceil(tmax/60)],[1 1],':k');             % dotted line at one
xlabel('Time (min)');
ylabel('Relative R');

% Green

figure;  hold on;
title(sprintf('Green display values vs time, wavelength = %d nm',filter_lambda));
GS = zeros(tmax,LN,'double');
gvi = zeros(tmax,LN,'double');

for n = 1:LN
  gv = squeeze(RGBmatch(2,n,:));       % G values for one lightness, all iterations
  gm = mean(gv);                       % mean over iterations
  th = squeeze(Tmatch(1,n,:));
  tm = squeeze(Tmatch(2,n,:));
  ts = squeeze(Tmatch(3,n,:));
  te = 3600*(th-t1h) + 60*(tm-t1m) + (ts-t1s);  % elapsed time for each observation (sec)
  GS(:,n) = interp1(te,gv,1:tmax,'linear','extrap');  % interpolate to seconds  
  gvi(:,n) = GS(:,n)'/gm;                   % normalise to average over all iterations
  plot(te/60,gv/gm,'-k');
end
gvm = mean(gvi(:,2:LN-1),2);                % mean over all lightnesses except first and last
plot(trange/60,gvm,'-g','LineWidth',3);     % plot as thick red line
plot([1 ceil(tmax/60)],[1 1],':k');         % dotted line at one
xlabel('Time (min)');
ylabel('Relative G');

% Blue

figure;  hold on;
title(sprintf('Blue display values vs time, wavelength = %d nm',filter_lambda));
BS = zeros(tmax,LN,'double');
bvi = zeros(tmax,LN,'double');

for n = 1:LN
  bv = squeeze(RGBmatch(3,n,:));       % B values for one lightness, all iterations
  bm = mean(bv);                       % mean over iterations
  th = squeeze(Tmatch(1,n,:));
  tm = squeeze(Tmatch(2,n,:));
  ts = squeeze(Tmatch(3,n,:));
  te = 3600*(th-t1h) + 60*(tm-t1m) + (ts-t1s);  % elapsed time for each observation (sec)
  BS(:,n) = interp1(te,bv,1:tmax,'linear','extrap');  % interpolate to seconds  
  bvi(:,n) = BS(:,n)'/bm;                   % normalise to average over all iterations
  plot(te/60,bv/bm,'-k');
end
bvm = mean(bvi(:,2:LN-1),2);                % mean over all lightnesses except first and last
plot(trange/60,bvm,'-b','LineWidth',3);     % plot as thick red line
plot([1 ceil(tmax/60)],[1 1],':k');         % dotted line at one
xlabel('Time (min)');
ylabel('Relative B');

%% Compute CCT of match vs time

% Read standard observer CMFs (380-780 nm, 1nm intervals)

ciedir = fullfile('D:','Research at UCL','Colour standards','CIE colorimetric data');
cmffile = fullfile(ciedir,'StdObs-2deg-1nm.txt');
format = '%d %f %f %f';
fid = fopen(cmffile, 'r');
[Ar,count] = fscanf(fid, format, [4, inf]);  % read the whole file into array A
fclose(fid);

s1count = 401;
Xcmf = Ar(2,1:s1count);           % use range 380-780 nm
Ycmf = Ar(3,1:s1count);
Zcmf = Ar(4,1:s1count);
clear Ar

figure;  hold on;
title('CIE tristimulus functions');
plot(380:780,Xcmf,'-r');
plot(380:780,Ycmf,'-g');
plot(380:780,Zcmf,'-b');

% Plot display spectra

figure;  hold on;
title('LCD display primaries');
xlabel('Wavelength (nm)');
ylabel('Relative power');
plot(lambda,RGBW_spectrum(:,1),'-r');
plot(lambda,RGBW_spectrum(:,2),'-g');
plot(lambda,RGBW_spectrum(:,3),'-b');
plot(lambda,RGBW_spectrum(:,4),':k');

% Interpolate spectra of display primaries to 1nm intervals

RGB_SPD = zeros(s1count,4,'double');

for k = 1:4
  RGB_SPD(:,k) = interp1(lambda,RGBW_spectrum(:,k),380:780,'spline');
end

% Make lookup tables for display tone curves

RGBlut = zeros(3,256,'double');

for k = 1:3
  Lum = squeeze(XYZ(2,:,k))/XYZ(2,21,k);
  RGBlut(k,:) = interp1(sval,Lum,0:255,'spline'); % interpolate luminance
end

figure;  hold on;
title('Normalised luminance for display R,G,B');
plot(0:255,RGBlut(1,:),'-r');
plot(0:255,RGBlut(2,:),'-g');
plot(0:255,RGBlut(3,:),'-b');

% Calculate tristimulus values of display white reference

DW = squeeze(RGB_SPD(:,4))';
Norm = 100/sum(DW.*Ycmf);              % normalising factor
Xw = sum(DW.*Xcmf)*Norm;
Yw = sum(DW.*Ycmf)*Norm;               % calculate white reference
Zw = sum(DW.*Zcmf)*Norm;
fprintf('Display white XYZ = %5.3f,%5.3f,%5.3f\n',Xw,Yw,Zw);

%% Construct spectrum for RGB colour match value

LN = 8;                              % selected lightness value

nmin = floor((tmax-1)/60);           % number of minutes in session
LD = zeros(nmin,1,'double');
AD = zeros(nmin,1,'double');
BD = zeros(nmin,1,'double');

for m = 1:nmin
  t = 60*(m-1)+1;                    % index to seconds
  r = uint8(255*RS(t,LN));           % get R,G,B digital values
  g = uint8(255*GS(t,LN));
  b = uint8(255*BS(t,LN));
  rgb = double([r g b]);
  SPD = zeros(s1count,1,'double');
  for k = 1:3
    SPD = SPD+squeeze(RGB_SPD(:,k))*RGBlut(k,rgb(k));  % make composite spectrum
  end
  X = sum(SPD'.*Xcmf)*Norm;
  Y = sum(SPD'.*Ycmf)*Norm;           % calculate tristimulus values
  Z = sum(SPD'.*Zcmf)*Norm;
  [LD(m) AD(m) BD(m)] = XYZtoLAB(X,Y,Z,Xw,Yw,Zw);  % convert to LAB
end

%figure;  hold on;
%title('Composite SPD');
%plot(380:780,SPD,'-k');

%% Plot a*,b* as a function of time

figure;  hold on;  grid on;  axis square;  rotate3d;
title(sprintf('a*b* values vs time, wavelength = %d nm',filter_lambda));
xlabel('a*');
ylabel('b*');
zlabel('Time (minutes)');
plot3(AD,BD,1:m,'-k','LineWidth',2);
for m = 1:nmin
  t = 60*(m-1)+1;                    % index to seconds
  r = RS(t,LN);                      % get R,G,B digital values
  g = GS(t,LN);
  b = BS(t,LN);
  rgb = double([r g b]);
  plot3(AD(m),BD(m),m,'ok','MarkerSize',10,'MarkerEdgeColor','none',...
      'MarkerFaceColor',rgb);
end
