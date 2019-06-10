%-------------------------------------------------------------------------
%  VISION SPHERE EXPERIMENT - WHITE BALANCE - TIME SERIES
%
%  Analyses experimental results by wavelength variation and plot figures
%
%-------------------------------------------------------------------------

dir = 'Lindsay time series - Sep 2012';      % data directory
%dir = 'Tania  time series - Apr 2013';      % data directory
wmin = 400;  wmax = 700;            % range of wavelengths (20nm intervals)
wrange = wmin:20:wmax;
N = 10;                             % number of repetitions over time
LN = 16;                            % number of lightness levels per repeat
WN = length(wrange);                % number of wavelength samples

LABM = zeros(3,LN,N,WN,'double');
RGBM = zeros(3,LN,N,WN,'double');
TM = zeros(3,LN,N,WN,'double');

for n = 1:WN
  lambda = wrange(n);               % filter wavelength
  fname = sprintf('%dnm - time',lambda);
  load(fullfile(dir,fname));            % load experimental results for wavelength
  LABM(:,:,:,n) = LABmatch(:,1:LN,:);
  RGBM(:,:,:,n) = RGBmatch(:,1:LN,:);   % copy all data into a single array
  TM(:,:,:,n) = Tmatch(:,1:LN,:);
end

Lval = squeeze(LABM(1,:,1,1));      % L values

%% Read display data for primaries and luminance calibation

fname = 'LCD display measurement.mat';
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

Rlut = squeeze(RGBlut(1,:));        % extract 1D versions
Glut = squeeze(RGBlut(2,:));
Blut = squeeze(RGBlut(3,:));

rw = RGBW_XYZ(2,1);                 % luminance of R,G,B display primaries (max signal)
gw = RGBW_XYZ(2,2);
bw = RGBW_XYZ(2,3);

% Convert match values to display luminance

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

% Plot three tone curves

figure;  hold on;
title('Display luminance vs input signal');
plot(ndac,RGBlut(1,:),'-r');
plot(ndac,RGBlut(2,:),'-g');
plot(ndac,RGBlut(3,:),'-b');
xlabel('Input signal (normalised)');
ylabel('Screen luminance (cd/m^2)');

%% Calculate session lengths for each wavelength

TI = zeros(4,WN,'double');          % start time (h,m,s) and elapsed time (sec)

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
   mv = mean(mean(vv));             % mean over all iterations and wavelengths
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

%% Plot all observed colours in RGB display values

p = 250;                            % position of axis markers
m = 10;                             % size of axis markers

figure;  hold on;  grid on;  rotate3d on;  axis equal;
title('RGB display values for all phases');

for w = 1:WN
  for n = 1:N
    for i = 2:LN-1
       rgb = RGBM(:,i,n,w);         % match values (range 0-1 if in gamut)
       RGB = rgb*255;               % rescale to 8-bit range
       if ((rgb(1)<0)||(rgb(2)<0)||(rgb(3)<0)) rgb = [1 0 0]; end
       if ((rgb(1)>1)||(rgb(2)>1)||(rgb(3)>1)) rgb = [1 0 0]; end
       plot3(RGB(1),RGB(2),RGB(3),'o','MarkerSize',5,'MarkerEdgeColor',rgb,...
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
plot3(p,0,0,'or','MarkerSize',m,'MarkerFaceColor','r');
plot3(0,p,0,'og','MarkerSize',m,'MarkerFaceColor','g');
plot3(0,0,p,'ob','MarkerSize',m,'MarkerFaceColor','b');
plot3(0,p,p,'oc','MarkerSize',m,'MarkerFaceColor','c');
plot3(p,0,p,'om','MarkerSize',m,'MarkerFaceColor','m');
plot3(p,p,0,'oy','MarkerSize',m,'MarkerFaceColor','y');
plot3(p,p,p,'ok','MarkerSize',m,'MarkerFaceColor','w');
xlabel('R');
ylabel('G');
zlabel('B');

%% Plot all observed colours in RGB luminance

p = 1;                              % position of axis markers
m = 10;                             % size of axis markers

figure;  hold on;  grid on;  rotate3d on;  axis equal;
title('RGB display luminance for all phases');

for w = 1:WN
  for n = 1:N
    for i = 1:LN
       rgb = RGBL(:,i,n,w);         % normalised luminance values (range 0-1)
       T = rgb>1;  rgb(T) = 1;      % clip values to range 0-1
       plot3(rgb(1),rgb(2),rgb(3),'o','MarkerSize',5,'MarkerEdgeColor',rgb,...
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
plot3(p,0,0,'or','MarkerSize',m,'MarkerFaceColor','r');
plot3(0,p,0,'og','MarkerSize',m,'MarkerFaceColor','g');
plot3(0,0,p,'ob','MarkerSize',m,'MarkerFaceColor','b');
plot3(0,p,p,'oc','MarkerSize',m,'MarkerFaceColor','c');
plot3(p,0,p,'om','MarkerSize',m,'MarkerFaceColor','m');
plot3(p,p,0,'oy','MarkerSize',m,'MarkerFaceColor','y');
plot3(p,p,p,'ok','MarkerSize',m,'MarkerFaceColor','w');
xlabel('Normalised R luminance');
ylabel('Normalised G luminance');
zlabel('Normalised B luminance');

%% Plot LAB slider settings for all observed colours

f = 50;                             % a,b axis scaling factor

figure;  hold on;  grid on;  rotate3d on;
title('LAB slider settings for all phases');

for w = 1:WN
  for n = 1:N
    for i = 1:LN
       ld = LABM(1,i,n,w);
       ad = LABM(2,i,n,w);
       bd = LABM(3,i,n,w);
       rgb = double(RGBD(:,i,n,w))/255;
       plot3(ad,bd,ld,'o','MarkerSize',5,'MarkerEdgeColor',rgb,...
      'MarkerFaceColor',rgb);
    end
  end
end

axis([-f f -f f 0 100]);
plot3([0 0],[-f f],[50 50],':k');
plot3([-f f],[0 0],[50 50],':k');   % draw axes
plot3([0 0],[0 0],[0 100],':k');
plot3(f-p,0,50,'or','MarkerSize',m,'MarkerFaceColor','r');
plot3(0,f-p,50,'oy','MarkerSize',m,'MarkerFaceColor','y');
plot3(-(f-p),0,50,'og','MarkerSize',m,'MarkerFaceColor','g');
plot3(0,-(f-p),50,'ob','MarkerSize',m,'MarkerFaceColor','b');
plot3(0,0,p,'ok','MarkerSize',m,'MarkerFaceColor','k');
plot3(0,0,50,'ok','MarkerSize',m,'MarkerFaceColor',[0.5 0.5 0.5]);
plot3(0,0,100-p,'ok','MarkerSize',m,'MarkerFaceColor','w');
xlabel('A');
ylabel('B');
zlabel('L');

%% Plot all observed colours in LAB coordinates under D65

f = 30;                             % a,b axis scaling factor
p = 2;                              % inset distance for axis markers
m = 8;                              % size of axis markers

figure;  hold on;  grid on;  rotate3d on;  axis equal;
title('CIELAB coordinates for all phases, converted from display RGB');

for w = 1:WN
  for n = 1:N
    for i = 1:LN-3
       ld = LABD(1,i,n,w);
       ad = LABD(2,i,n,w);
       bd = LABD(3,i,n,w);
       rgb = double(RGBD(:,i,n,w))/255;
       plot3(ad,bd,ld,'o','MarkerSize',5,'MarkerEdgeColor',rgb,...
      'MarkerFaceColor',rgb);
    end
  end
end

axis([-f f -f+10 f+10 20 80]);
plot3([0 0],[-f f],[50 50],':k');
plot3([-f f],[0 0],[50 50],':k');   % draw axes
plot3([0 0],[0 0],[0 100],':k');
plot3(f-p,0,50,'or','MarkerSize',m,'MarkerFaceColor','r');
plot3(0,f+10-p,50,'oy','MarkerSize',m,'MarkerFaceColor','y');
plot3(-f+p,0,50,'og','MarkerSize',m,'MarkerFaceColor','g');
plot3(0,-f+10+p,50,'ob','MarkerSize',m,'MarkerFaceColor','b');
plot3(0,0,20+p,'ok','MarkerSize',m,'MarkerFaceColor','k');
plot3(0,0,50,'ok','MarkerSize',m,'MarkerFaceColor',[0.5 0.5 0.5]);
plot3(0,0,80-p,'ok','MarkerSize',m,'MarkerFaceColor','w');
xlabel('a*');
ylabel('b*');
zlabel('L*');

%% RGB vs wavelength for fixed values of iteration and lightness

iter = 6;
lness = 10;

figure;  hold on;
title('Raw RGB values for neutral match');
L = Lval(lness);
rmax = max(RGBM(1,lness,:));
text(wmin+20,rmax-0.03,sprintf('Iteration = %d, Lightness = %d',iter,L));
rs = squeeze(RGBM(1,lness,iter,:));
plot(wrange,rs,'-r');                   % R
gs = squeeze(RGBM(2,lness,iter,:));
plot(wrange,gs,'-g');                   % G
bs = squeeze(RGBM(3,lness,iter,:));
plot(wrange,bs,'-b');                   % B
xlabel('Filter wavelength (nm)');
ylabel('Normalised R,G,B (range 0-1)');

%% 3D surface vs wavelength and iteration

lness = 8;
[XL,YL] = meshgrid(1:N,wrange);       % iteration on X axis; wavelength on Y axis

figure;  hold on;  rotate3d;  grid on;
RL = squeeze(RGBM(1,lness,:,:));      % red display values for fixed lightness
ZL = Rlut(round(imax*RL)+1);          % convert to luminance values
surf(XL,YL,ZL');
title(sprintf('Lightness %d',Lval(lness)));
xlabel('Iteration');
ylabel('Wavelength (nm)');
zlabel('Red display luminance (cd/m^2)');

figure;  hold on;  rotate3d;  grid on;
GL = squeeze(RGBM(2,lness,:,:));      % green display values for fixed lightness
ZL = Glut(round(imax*GL)+1);          % convert to luminance values
surf(XL,YL,ZL');
title(sprintf('Lightness %d',Lval(lness)));
xlabel('Iteration');
ylabel('Wavelength (nm)');
zlabel('Green display luminance (cd/m^2)');

figure;  hold on;  rotate3d;  grid on;
BL = squeeze(RGBM(3,lness,:,:));      % blue display values for fixed lightness
ZL = Blut(round(imax*BL)+1);          % convert to luminance values
surf(XL,YL,ZL');
title(sprintf('Lightness %d',Lval(lness)));
xlabel('Iteration');
ylabel('Wavelength (nm)');
zlabel('Blue display luminance (cd/m^2)');

%% 3D surface vs wavelength and lightness

iter = 5;
[XL,YL] = meshgrid(Lval(1:LN),wrange);  % lightness on X axis; wavelength on Y axis

figure;  hold on;  rotate3d;  grid on;
title(sprintf('Display RGB vs wavelength and lightness, iteration = %d',iter));
RL = squeeze(RGBM(1,:,iter,:));       % red display values for fixed iteration
T = RL<0;  RL(T) = 0;                 % clip to minimum
T = RL>1;  RL(T) = 1;                 % clip to maximum
ZL = Rlut(round((imax-1)*RL)+1);      % convert to luminance values
ZM = mean(ZL,2);                      % mean over all wavelengths at each lightness
for i = 1:LN
  ZL(i,:) = ZL(i,:)/ZM(i);            % normalise to mean lightness over all wavelengths
end
surf(XL,YL,ZL');
axis([20 85 400 700 min(min(ZL)) max(max(ZL))]);
xlabel('Lightness');
ylabel('Wavelength (nm)');
zlabel('Relative red display luminance');

figure;  hold on;  rotate3d;  grid on;
title(sprintf('Display RGB vs wavelength and lightness, iteration = %d',iter));
GL = squeeze(RGBM(2,:,iter,:));       % green display values for fixed iteration
T = GL<0;  GL(T) = 0;                 % clip to minimum
T = GL>1;  GL(T) = 1;                 % clip to maximum
ZL = Glut(round((imax-1)*GL)+1);      % convert to luminance values
ZM = mean(ZL,2);                      % mean over all wavelengths at each lightness
for i = 1:LN
  ZL(i,:) = ZL(i,:)/ZM(i);            % normalise to mean lightness over all wavelengths
end
surf(XL,YL,ZL');
axis([20 85 400 700 min(min(ZL)) max(max(ZL))]);
xlabel('Lightness');
ylabel('Wavelength (nm)');
zlabel('Relative green display luminance');

figure;  hold on;  rotate3d;  grid on;
title(sprintf('Display RGB vs wavelength and lightness, iteration = %d',iter));
BL = squeeze(RGBM(3,:,iter,:));       % blue display values for fixed iteration
T = BL<0;  BL(T) = 0;                 % clip to minimum
T = BL>1;  BL(T) = 1;                 % clip to maximum
ZL = Blut(round((imax-1)*BL)+1);      % convert to luminance values
ZM = mean(ZL,2);                      % mean over all wavelengths at each lightness
for i = 1:LN
  ZL(i,:) = ZL(i,:)/ZM(i);            % normalise to mean lightness over all wavelengths
end
surf(XL,YL,ZL');
axis([20 85 400 700 min(min(ZL)) max(max(ZL))]);
xlabel('Lightness');
ylabel('Wavelength (nm)');
zlabel('Relative blue display luminance');

%% Plot RGB for one wavelength

lam = 6;
lambda = wrange(lam);

figure;  hold on;
title(sprintf('Display RGB vs time, adapting wavelength = %d nm',lambda));
tmax = TI(4,lam);                            % session length (sec) for this wavelength
xlim([0 tmax/60]);                           % X axis in minutes
rvi = squeeze(RGBrel(1,:,1:tmax,lam));       % extract all red values for wavelength
rvm = mean(rvi(1:LN,:));                     % mean over all lightnesses
plot((1:tmax)/60,rvm,'-r','LineWidth',2);    % plot as thick line
gvi = squeeze(RGBrel(2,:,1:tmax,lam));       % extract all green values for wavelength
gvm = mean(gvi(1:LN,:));                     % mean over all lightnesses
plot((1:tmax)/60,gvm,'-g','LineWidth',2);    % plot as thick line
bvi = squeeze(RGBrel(3,:,1:tmax,lam));       % extract all blue values for wavelength
bvm = mean(bvi(1:LN,:));                     % mean over all lightnesses
plot((1:tmax)/60,bvm,'-b','LineWidth',2);    % plot as thick line
tm = ceil(max(TI(4,:))/60);                  % length of longest session (minutes)
plot([1 tm],[1 1],':k');                     % dotted line at one
xlabel('Time (min)');
ylabel('Value relative to mean over all iterations');

%% Plot surface of R/B ratio for wavelength vs time

lness = 8;
ti = 1;                 % time interval (mins)
tmax = floor(TNM/ti);   % number of time steps
wi = 400:5:700;
wl = length(wi);

RBrat = zeros(wl,tmax,'double');

for tm = 1:tmax
  t = tm*ti;                              % time (mins)
  r = squeeze(RGBval(1,lness,t*60,:));    % extract red values for all wavelengths
  ri = interp1(wrange,r,wi,'spline');
  b = squeeze(RGBval(3,lness,t*60,:));    % extract blue values for all wavelengths
  bi = interp1(wrange,b,wi,'spline');
  RBrat(:,tm) = ri./bi;
end

[XL,YL] = meshgrid(wi,1:tmax);

figure;  hold on;  rotate3d;  grid on;
title(sprintf('Display R/B ratio for wavelength vs time, lightness = %d',Lval(lness)));
ylim([0 TNM]);
surf(XL,YL,RBrat');
xlabel('Wavelength (nm)');
ylabel('Time (min)');
zlabel('R/B ratio');

%% Image mosaic of matching colours vs wavelength and time (5 minute intervals)

tint = 5;                       % time interval (min)
TN = floor(TNM/tint);           % number of time steps
b = 40;                         % pixels in box side
s = 4;                          % spacing between boxes
w = s+WN*(b+s);                 % width of array (wavelength axis)
h = s+LN*(b+s);                 % height of array (lightness axis)
Im = zeros(h,w,3,'uint8');      % image array
idir = fullfile('C:','Research at UCL','Matlab scripts','Experiment','Figures',...
    'Colour mosaics by time period');

for tm = 1:TN
 t = tint*tm;                   % time in minutes
 for n = 1:WN
  xp = s+(n-1)*(b+s);           % x pixel address
  for i = 1:LN
    rgb = RGBD(:,i,t,n);        % get sRGB value (display signal, 8-bit)
    yp = s+(i-1)*(b+s);         % y pixel address
    Im(yp:yp+b-1,xp:xp+b-1,1) = rgb(1);  % fill one square in array
    Im(yp:yp+b-1,xp:xp+b-1,2) = rgb(2);
    Im(yp:yp+b-1,xp:xp+b-1,3) = rgb(3);
  end
 end
 iname = fullfile(idir,sprintf('%d_mins.tif',t));
% imwrite(Im,iname,'tif');           % write the image
end

%% Image mosaic of matching colours vs wavelength and lightness

w = s+WN*(b+s);                 % width of array (wavelength axis)
h = s+TN*(b+s);                 % height of array (time axis)
Im = zeros(h,w,3,'uint8');      % image array
idir = fullfile('C:','Research at UCL','Matlab scripts','Experiment','Figures',...
    'Colour mosaics by lightness');

 for i = 1:LN                   % separate grid for each lightness
  for n = 1:WN                  % wavelength on x axis
   xp = s+(n-1)*(b+s);          % x pixel address
   for tm = 1:TN                % time interval on y axis
    t = tm*tint;
    rgb = RGBD(:,i,t,n);        % get sRGB value (display signal, 8-bit)
    yp = s+(tm-1)*(b+s);         % y pixel address
    Im(yp:yp+b-1,xp:xp+b-1,1) = rgb(1);  % fill one square in array
    Im(yp:yp+b-1,xp:xp+b-1,2) = rgb(2);
    Im(yp:yp+b-1,xp:xp+b-1,3) = rgb(3);
  end
 end
 iname = fullfile(idir,sprintf('Lightness_%d.tif',Lval(i)));
% imwrite(Im,iname,'tif');           % write the image
end

%% Plot as 3D surface

Csurf = zeros(TN,WN,'double');
Cmap = zeros(TN,WN,3,'double');
[XL,YL] = meshgrid(wrange,tint:tint:tint*TN);

i = 5;          % lightness (constant)

for n = 1:WN
  for t = 1:TN
    rgb = RGBD(:,i,t*tint,n);   % get RGB value (display signal, 8-bit)
    Cmap(t,n,:) = double(rgb)/255;
  end
end

figure;  hold on;  grid on;  rotate3d;  axis square;
title(sprintf('Lightness %d',i));
Csurf(:,:) = 1;
surf(XL,YL,Csurf,Cmap);
xlabel('Wavelength (nm)');
ylabel('Time elapsed (min)');
zlabel('');
