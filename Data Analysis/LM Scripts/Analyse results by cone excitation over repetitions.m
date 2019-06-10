%-------------------------------------------------------------------------
%  VISION SPHERE EXPERIMENT - WHITE BALANCE - TIME SERIES
%
%  Analyse experimental results in terms of L,M,S cone excitations
%  Compare results over multiple repetitions at same wavelength
%
%-------------------------------------------------------------------------

dir = fullfile('C:','Research at UCL','Experiment','Tania  time series - Apr 2013');
s1 = 341;               % number of values in 1nm spectrum 390-730 nm
w1 = 390:1:730;
s5 = 69;                % number of values in 5nm spectrum 390-730 nm
w5 = 390:5:730;
ciedir = fullfile('C:','Research at UCL','Colour standards','CIE colorimetric data');
cvrldir = fullfile('C:','Research at UCL','Colour standards','CVRL cone fundamentals');

%% Read Stockman-Sharpe 10-deg cone fundamentals (390-730 nm in 1nm intervals)

ssfile = fullfile(cvrldir,'Stockman-Sharpe cone fundamentals - lin-2deg-1nm.txt');
format = '%d %f %f %f';
fid = fopen(ssfile,'r');
[Obs,count] = fscanf(fid,format,[4,inf]);  % read the whole file into array
fclose(fid);

Lcone = Obs(2,1:s1);               % extract cone response data fields
Mcone = Obs(3,1:s1);               % data range 390-730nm
Scone = Obs(4,1:s1);
LMScone = [Lcone; Mcone; Scone];        % 3x341 array

%% Plot cone fundamentals

figure;  hold on;
title('CVRL 2-deg cone fundamentals');
plot(w1,Lcone,'-r');
plot(w1,Mcone,'-g');
plot(w1,Scone,'-b');
legend('L cone','M cone','S cone');
xlabel('Wavelength (nm)');
ylabel('Relative sensitivity');

%% Read experimental data for neutral matching

lambda = 500;                       % filter wavelength
N = 10;                             % number of iterations over time
LN = 16;                            % number of lightness levels per repeat
RN = 5;                             % number of iterations at same wavelength
LABM = zeros(3,LN,N,RN,'double');
RGBM = zeros(3,LN,N,RN,'double');
TM = zeros(3,LN,N,RN,'double');

for n = 1:RN
  fname = fullfile(dir,sprintf('%dnm - time',lambda));
  load(fname);                      % load experimental results for wavelength
  LABM(:,:,:,n) = LABmatch(:,1:LN,:);
  RGBM(:,:,:,n) = RGBmatch(:,1:LN,:);   % copy all data into a single array
  TM(:,:,:,n) = Tmatch(:,1:LN,:);
end

Lval = squeeze(LABM(1,:,1,1));      % L values

%% Read display data for primaries and luminance calibation

dfile = fullfile('C:','Research at UCL','Experiment','Large LCD display measurement.mat');
load(dfile);                        % load display data

imax = 65535;                       % maximum integer in 16-bit range
dac = 0:imax;                       % range of input values
ndac = dac/imax;                    % normalised to range 0-1
RGBlut = zeros(3,imax+1,'double');  % lookup tables for tone curve (quasi 16-bit)
RGBL = zeros(3,LN,N,RN,'double');   % converted normalised luminance values

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

for n = 1:RN
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
  DS(:,k) = interp1(w4,p,w1','spline');  % interpolate to range 390-730nm
end

figure;  hold on;
title('Display RGB spectra');
plot(w1,DS(:,1),'-r');
plot(w1,DS(:,2),'-g');
plot(w1,DS(:,3),'-b');
xlabel('Wavelength (nm)');
ylabel('Power');

%% Calculate LMS response to RGB display primaries

RGB_LMS = zeros(3,3,'double');  % RGB columns, LMS rows

for i = 1:3
  for k = 1:3
    RGB_LMS(k,i) = sum(DS(:,i)'.*LMScone(k,:));  % sum over range 390-730nm
  end
end

RGB_LMS = RGB_LMS/RGB_LMS(2,2);     % normalise to green on M cone

%% Analyse relative RGB vs time

TI = zeros(4,RN,'double');          % start time (h,m,s) and elapsed time (sec)

% Calculate session lengths for each wavelength

for n = 1:RN
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
RGBval = zeros(3,LN,TNS,RN,'double');      % actual value per second (range 0-1)
RGBrel = zeros(3,LN,TNS,RN,'double');      % value relative to mean over repetitions

for n = 1:RN
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
XYZD = zeros(3,LN,TNM,RN,'double');     % display XYZ of match
XYZD65 = zeros(3,LN,TNM,RN,'double');   % display XYZ of match (D65 white ref)
LABD = zeros(3,LN,TNM,RN,'double');     % display LAB of match
RGBD = zeros(3,LN,TNM,RN,'uint8');      % display RGB of match (sRGB, D65)

for n = 1:RN
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

%% Convert display RGB to LMS cone excitations at 1 minute intervals

rgbw = RGBW_XYZ(2,:);     % luminance of R,G,B display primaries (max signal)

TNM = floor(min(TI(4,:)/60));       % find length of shortest session (minutes)
LMS = zeros(3,LN,TNM,RN,'double');  % LMS of match

for n = 1:RN
 for t = 1:TNM
  for i = 1:LN
    ts = 60*t;                  % time in seconds
    rgb = RGBval(:,i,ts,n);     % get RGB match value (display signal, range 0-1)
    rgbi = 1+uint16(65534*rgb); % convert to 16-bit LUT index
    rl = Rlut(rgbi(1))/rw;      % convert to fraction of white luminance
    gl = Glut(rgbi(2))/gw;
    bl = Blut(rgbi(3))/bw;
    rgbl = [rl gl bl];          % RGB luminances
    for k = 1:3
      lmsvec = RGB_LMS(k,:)';   % row of matrix for one of L,M,S
      LMS(k,i,t,n) = rgbl*lmsvec;   % sum contribution from each primary
    end
  end
 end
end

%% Plot cone excitations in LMS cube

p = 1;                              % position of axis markers
m = 10;                             % size of axis markers

lmsmax = max(max(max(max(LMS))));   % max value for scaling

figure;  hold on;  grid on;  rotate3d on;  axis equal;
title('LMS cone excitation values for all matches');

for w = 1:RN
  for n = 1:N
    for i = 1:LN
       lms = LMS(:,i,n,w)/lmsmax;        % normalised cone values (range 0-1)
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

%% Plot LMS vs repetition and time for constant lightness

lness = 8;              % lightness (fixed)
RI = 1:5;

lp = reshape(LMS(1,lness,:,:),TNM*RN,1);    % extract L as vector
ml = mean(lp);                              % mean L over plane
sl = std(lp);
lm = reshape(LMS(2,lness,:,:),TNM*RN,1);    % extract M as vector
mm = mean(lm);                              % mean M over plane
sm = std(lm);
ls = reshape(LMS(3,lness,:,:),TNM*RN,1);    % extract S as vector
ms = mean(ls);                              % mean S over plane
ss = std(ls);
fprintf(1,'\nStd and mean cone excitations for target lightness %d\n',Lval(lness));
fprintf(1,'L %f/%f=%5.3f,  M %f/%f=%5.3f,  S %f/%f=%5.3f\n',...
            sl,ml,sl/ml,sm,mm,sm/mm,ss,ms,ss/ms); 

% Plot figures

figure;  hold on;  grid on;
title(sprintf('L cone vs time for five repetitions, lightness = %d',Lval(lness)));
for n = 1:RN
  zl = squeeze(LMS(1,lness,:,n));       % extract all values for one run
  plot(1:TNM,zl,'-r');
end
xlabel('Time (min)');
ylabel('L cone excitation');
%axis([1 TNM 0.35 0.45]);

figure;  hold on;  grid on;
title(sprintf('M cone vs time for five repetitions, lightness = %d',Lval(lness)));
for n = 1:RN
  zl = squeeze(LMS(2,lness,:,n));       % extract all values for one run
  plot(1:TNM,zl,'-g');
end
xlabel('Time (min)');
ylabel('M cone excitation');
%axis([1 TNM 0.3 0.4]);

figure;  hold on;  grid on;
title(sprintf('S cone vs time for five repetitions, lightness = %d',Lval(lness)));
for n = 1:RN
  zl = squeeze(LMS(3,lness,:,n));       % extract all values for one run
  plot(1:TNM,zl,'-b');
end
xlabel('Time (min)');
ylabel('S cone excitation');
%axis([1 TNM 0.15 0.35]);
