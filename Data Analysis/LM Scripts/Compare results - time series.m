%--------------------------------------------------------------------------
%  VISION SPHERE EXPERIMENT - WHITE BALANCE - TIME SERIES
%
%  Compare two sets of experimental results for one wavelength
%
%--------------------------------------------------------------------------

filter_lambda = 500;                       % filter wavelength
dir1 = 'Results - Apr 2013';
name1 = 'TR1';
dir2 = 'Results - Apr 2013';
name2 = 'TR2';
fname1 = sprintf('Tania %dnm - time.mat',filter_lambda);
fname2 = sprintf('Tania %dnm - time v2.mat',filter_lambda);

%load(fullfile(dir1,fname1));        % load experimental results for wavelength
load(fname1);
LABmatch1 = LABmatch;
RGBmatch1 = RGBmatch;
Tmatch1 = Tmatch;
Wref1 = Wref;
fprintf('\nDataset #1 = %s\n',name1);

%load(fullfile(dir2,fname2));
load(fname2);
LABmatch2 = LABmatch;
RGBmatch2 = RGBmatch;
Tmatch2 = Tmatch;
Wref2 = Wref;
fprintf('Dataset #2 = %s\n',name2);

if Wref1(1) ~= Wref2(1)
  disp('ERROR - different display white reference');
end

N = 10;                             % number of repetitions over time
LN = 16;                            % number of lightness levels per repeat

Lval = squeeze(LABmatch1(1,:,1));    % L values (same for both datasets)

load('Large LCD display measurement.mat');  % load display data

%% Make mosaics of matching colours - lightness vs iteration

b = 40;                         % pixels in box side
s = 4;                          % spacing between boxes
w = s+N*(b+s);                  % width of array (iteration axis)
h = s+LN*(b+s);                 % height of array (lightness axis)

Im1 = zeros(h,w,3,'uint8');      % image array
for n = 1:N
  xp = s+(n-1)*(b+s);                    % x pixel address (iteration axis)
  for i = 1:LN
    rs = RGBmatch1(1,i,n);                % get R value (display signal, 8-bit)
    gs = RGBmatch1(2,i,n);
    bs = RGBmatch1(3,i,n);
    yp = s+(i-1)*(b+s);                  % y pixel address (lightness axis)
    Im1(yp:yp+b-1,xp:xp+b-1,1) = uint8(255*rs);  % fill one square in array
    Im1(yp:yp+b-1,xp:xp+b-1,2) = uint8(255*gs);
    Im1(yp:yp+b-1,xp:xp+b-1,3) = uint8(255*bs);
  end
end

figure;  hold on;
title(sprintf('%s data - Selected colours, wavelength = %d nm',name1,filter_lambda));
imshow(Im1);
xlabel('Iteration');
ylabel('Lightness');

Im2 = zeros(h,w,3,'uint8');      % image array
for n = 1:N
  xp = s+(n-1)*(b+s);                    % x pixel address (iteration axis)
  for i = 1:LN
    rs = RGBmatch2(1,i,n);                % get R value (display signal, 8-bit)
    gs = RGBmatch2(2,i,n);
    bs = RGBmatch2(3,i,n);
    yp = s+(i-1)*(b+s);                  % y pixel address (lightness axis)
    Im2(yp:yp+b-1,xp:xp+b-1,1) = uint8(255*rs);  % fill one square in array
    Im2(yp:yp+b-1,xp:xp+b-1,2) = uint8(255*gs);
    Im2(yp:yp+b-1,xp:xp+b-1,3) = uint8(255*bs);
  end
end

figure;  hold on;
title(sprintf('%s data - Selected colours, wavelength = %d nm',name2,filter_lambda));
imshow(Im2);
xlabel('Iteration');
ylabel('Lightness');

%% R/B ratio vs iteration

figure;  hold on;
title(sprintf('Ratio R/B display values vs iteration, wavelength = %d nm',filter_lambda));
plot([1 N],[1 1],':k');             % dotted line at one
xlabel('Iteration');
ylabel('R/B');

for n = 1:LN
  rv1 = squeeze(RGBmatch1(1,n,:));
  bv1 = squeeze(RGBmatch1(3,n,:));
  plot(1:N,rv1./bv1,':k');          % plot for individual iterations
end

for n = 1:N
  rv1 = squeeze(RGBmatch1(1,1:LN,n))';
  bv1 = squeeze(RGBmatch1(3,1:LN,n))';
  rbm1(n) = mean(rv1./bv1);
end
h1 = plot(1:N,rbm1,':m','LineWidth',3);   % plot mean of all iterations

for n = 1:LN
  rv2 = squeeze(RGBmatch2(1,n,:));
  bv2 = squeeze(RGBmatch2(3,n,:));
  plot(1:N,rv2./bv2,'-k');          % plot for individual iterations
end

for n = 1:N
  rv2 = squeeze(RGBmatch2(1,1:LN,n))';
  bv2 = squeeze(RGBmatch2(3,1:LN,n))';
  rbm2(n) = mean(rv2./bv2);
end
h2 = plot(1:N,rbm2,'-m','LineWidth',3);   % plot mean of all iterations

legend([h1 h2],name1,name2,'Location','NorthEast');

%% Interpolate relative R,G,B values vs time

t1h1 = squeeze(Tmatch1(1,1,1));               % time of first observation
t1m1 = squeeze(Tmatch1(2,1,1));
t1s1 = squeeze(Tmatch1(3,1,1));
t2h1 = squeeze(Tmatch1(1,LN,N));              % time of last observation
t2m1 = squeeze(Tmatch1(2,LN,N));
t2s1 = squeeze(Tmatch1(3,LN,N));
tmax1 = round(3600*(t2h1-t1h1)+60*(t2m1-t1m1)+(t2s1-t1s1));   % length of session from first to last (sec)
trange1 = 1:tmax1;
fprintf('\nSession length #1 = %d minutes\n',round(tmax1/60));
fprintf('Mean time per observation = %d seconds\n',round(tmax1/(LN*N)));

t1h2 = squeeze(Tmatch2(1,1,1));               % time of first observation
t1m2 = squeeze(Tmatch2(2,1,1));
t1s2 = squeeze(Tmatch2(3,1,1));
t2h2 = squeeze(Tmatch2(1,LN,N));              % time of last observation
t2m2 = squeeze(Tmatch2(2,LN,N));
t2s2 = squeeze(Tmatch2(3,LN,N));
tmax2 = round(3600*(t2h2-t1h2)+60*(t2m2-t1m2)+(t2s2-t1s2));   % length of session from first to last (sec)
trange2 = 1:tmax2;
fprintf('\nSession length #2 = %d minutes\n',round(tmax2/60));
fprintf('Mean time per observation = %d seconds\n',round(tmax2/(LN*N)));

% Red

figure;  hold on;
title(sprintf('Red display values vs time, wavelength = %d nm',filter_lambda));
plot([1 ceil(max(tmax1,tmax2)/60)],[1 1],':k');             % dotted line at one
xlabel('Time (min)');
ylabel('Relative R');

RS1 = zeros(tmax1,LN,'double');
rvi1 = zeros(tmax1,LN,'double');
for n = 1:LN
  rv = squeeze(RGBmatch1(1,n,:));       % R values for one lightness, all iterations
  rm = mean(rv);                       % mean over iterations
  th = squeeze(Tmatch1(1,n,:));
  tm = squeeze(Tmatch1(2,n,:));
  ts = squeeze(Tmatch1(3,n,:));
  te = 3600*(th-t1h1) + 60*(tm-t1m1) + (ts-t1s1);  % elapsed time for each observation (sec)
  RS1(:,n) = interp1(te,rv,1:tmax1,'linear','extrap');  % interpolate to seconds  
  rvi1(:,n) = RS1(:,n)'/rm;                   % normalise to average over all iterations
  plot(te/60,rv/rm,':k');
end
rvm1 = mean(rvi1(:,2:LN-1),2);                % mean over all lightnesses except first and last
h1 = plot(trange1/60,rvm1,':r','LineWidth',3);     % plot as thick red line
RS1(RS1<0) = 0;
RS1(RS1>1) = 1;

RS2 = zeros(tmax2,LN,'double');
rvi2 = zeros(tmax2,LN,'double');
for n = 1:LN
  rv = squeeze(RGBmatch2(1,n,:));       % R values for one lightness, all iterations
  rm = mean(rv);                       % mean over iterations
  th = squeeze(Tmatch2(1,n,:));
  tm = squeeze(Tmatch2(2,n,:));
  ts = squeeze(Tmatch2(3,n,:));
  te = 3600*(th-t1h2) + 60*(tm-t1m2) + (ts-t1s2);  % elapsed time for each observation (sec)
  RS2(:,n) = interp1(te,rv,1:tmax2,'linear','extrap');  % interpolate to seconds  
  rvi2(:,n) = RS2(:,n)'/rm;                   % normalise to average over all iterations
  plot(te/60,rv/rm,'-k');
end
rvm2 = mean(rvi2(:,2:LN-1),2);                % mean over all lightnesses except first and last
h2 = plot(trange2/60,rvm2,'-r','LineWidth',3);     % plot as thick red line
legend([h1 h2],name1,name2,'Location','NorthEast');
RS2(RS2<0) = 0;
RS2(RS2>1) = 1;

% Green

figure;  hold on;
title(sprintf('Green display values vs time, wavelength = %d nm',filter_lambda));
plot([1 ceil(max(tmax1,tmax2)/60)],[1 1],':k');             % dotted line at one
xlabel('Time (min)');
ylabel('Relative R');

GS1 = zeros(tmax1,LN,'double');
gvi1 = zeros(tmax1,LN,'double');
for n = 1:LN
  rv = squeeze(RGBmatch1(2,n,:));       % R values for one lightness, all iterations
  rm = mean(rv);                       % mean over iterations
  th = squeeze(Tmatch1(1,n,:));
  tm = squeeze(Tmatch1(2,n,:));
  ts = squeeze(Tmatch1(3,n,:));
  te = 3600*(th-t1h1) + 60*(tm-t1m1) + (ts-t1s1);  % elapsed time for each observation (sec)
  GS1(:,n) = interp1(te,rv,1:tmax1,'linear','extrap');  % interpolate to seconds  
  gvi1(:,n) = GS1(:,n)'/rm;                   % normalise to average over all iterations
  plot(te/60,rv/rm,':k');
end
gvm1 = mean(gvi1(:,2:LN-1),2);                % mean over all lightnesses except first and last
h1 = plot(trange1/60,gvm1,':g','LineWidth',3);     % plot as thick green line
GS1(GS1<0) = 0;
GS1(GS1>1) = 1;

GS2 = zeros(tmax2,LN,'double');
gvi2 = zeros(tmax2,LN,'double');
for n = 1:LN
  rv = squeeze(RGBmatch2(2,n,:));       % R values for one lightness, all iterations
  rm = mean(rv);                       % mean over iterations
  th = squeeze(Tmatch2(1,n,:));
  tm = squeeze(Tmatch2(2,n,:));
  ts = squeeze(Tmatch2(3,n,:));
  te = 3600*(th-t1h2) + 60*(tm-t1m2) + (ts-t1s2);  % elapsed time for each observation (sec)
  GS2(:,n) = interp1(te,rv,1:tmax2,'linear','extrap');  % interpolate to seconds  
  gvi2(:,n) = GS2(:,n)'/rm;                   % normalise to average over all iterations
  plot(te/60,rv/rm,'-k');
end
gvm2 = mean(gvi2(:,2:LN-1),2);                % mean over all lightnesses except first and last
h2 = plot(trange2/60,gvm2,'-g','LineWidth',3);     % plot as thick green line
legend([h1 h2],name1,name2,'Location','SouthEast');
GS2(GS2<0) = 0;
GS2(GS2>1) = 1;

% Blue

figure;  hold on;
title(sprintf('Blue display values vs time, wavelength = %d nm',filter_lambda));
plot([1 ceil(max(tmax1,tmax2)/60)],[1 1],':k');             % dotted line at one
xlabel('Time (min)');
ylabel('Relative R');

BS1 = zeros(tmax1,LN,'double');
bvi1 = zeros(tmax1,LN,'double');
for n = 1:LN
  rv = squeeze(RGBmatch1(3,n,:));       % R values for one lightness, all iterations
  rm = mean(rv);                       % mean over iterations
  th = squeeze(Tmatch1(1,n,:));
  tm = squeeze(Tmatch1(2,n,:));
  ts = squeeze(Tmatch1(3,n,:));
  te = 3600*(th-t1h1) + 60*(tm-t1m1) + (ts-t1s1);  % elapsed time for each observation (sec)
  BS1(:,n) = interp1(te,rv,1:tmax1,'linear','extrap');  % interpolate to seconds  
  bvi1(:,n) = BS1(:,n)'/rm;                   % normalise to average over all iterations
  plot(te/60,rv/rm,':k');
end
bvm1 = mean(bvi1(:,2:LN-1),2);                % mean over all lightnesses except first and last
h1 = plot(trange1/60,bvm1,':b','LineWidth',3);     % plot as thick blue line
BS1(BS1<0) = 0;
BS1(BS1>1) = 1;

BS2 = zeros(tmax2,LN,'double');
bvi2 = zeros(tmax2,LN,'double');
for n = 1:LN
  rv = squeeze(RGBmatch2(3,n,:));       % R values for one lightness, all iterations
  rm = mean(rv);                        % mean over iterations
  th = squeeze(Tmatch2(1,n,:));
  tm = squeeze(Tmatch2(2,n,:));
  ts = squeeze(Tmatch2(3,n,:));
  te = 3600*(th-t1h2) + 60*(tm-t1m2) + (ts-t1s2);  % elapsed time for each observation (sec)
  BS2(:,n) = interp1(te,rv,1:tmax2,'linear','extrap');  % interpolate to seconds  
  bvi2(:,n) = BS2(:,n)'/rm;                   % normalise to average over all iterations
  plot(te/60,rv/rm,'-k');
end
bvm2 = mean(bvi2(:,2:LN-1),2);                % mean over all lightnesses except first and last
h2 = plot(trange2/60,bvm2,'-b','LineWidth',3);     % plot as thick blue line
legend([h1 h2],name1,name2,'Location','NorthEast');
BS2(BS2<0) = 0;
BS2(BS2>1) = 1;

%% Compute LAB of match from reconstructed spectra

% Read standard observer CMFs (380-780 nm, 1nm intervals)

ciedir = fullfile('C:','Research at UCL','Colour standards','CIE colorimetric data');
cmffile = fullfile(ciedir,'StdObs-2deg-1nm.txt');
format = '%d %f %f %f';
fid = fopen(cmffile, 'r');
[Ar,count] = fscanf(fid, format, [4, inf]);  % read the whole file into array A
fclose(fid);

s1count = 401;
Xcmf = Ar(2,1:s1count);           % use wavelength range 380-780 nm
Ycmf = Ar(3,1:s1count);
Zcmf = Ar(4,1:s1count);
clear Ar

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

% Calculate tristimulus values of display white reference

DW = squeeze(RGB_SPD(:,4))';
Norm = 100/sum(DW.*Ycmf);              % normalising factor
Xw = sum(DW.*Xcmf)*Norm;
Yw = sum(DW.*Ycmf)*Norm;               % calculate white reference
Zw = sum(DW.*Zcmf)*Norm;
fprintf('Display white XYZ = %5.3f,%5.3f,%5.3f\n',Xw,Yw,Zw);

% Construct spectrum for RGB colour match value

nmin1 = floor((tmax1-1)/60);         % number of minutes in session
LD1 = zeros(nmin1,LN,'double');
AD1 = zeros(nmin1,LN,'double');
BD1 = zeros(nmin1,LN,'double');

for Lness = 1:LN
 for m = 1:nmin1
  t = 60*(m-1)+1;                    % index to seconds
  r = uint8(255*RS1(t,Lness));       % get R,G,B digital values
  g = uint8(255*GS1(t,Lness));
  b = uint8(255*BS1(t,Lness));
  rgb = double([r g b]);
  SPD = zeros(s1count,1,'double');
  for k = 1:3
    SPD = SPD+squeeze(RGB_SPD(:,k))*RGBlut(k,rgb(k));  % make composite spectrum
  end
  X = sum(SPD'.*Xcmf)*Norm;
  Y = sum(SPD'.*Ycmf)*Norm;           % calculate tristimulus values
  Z = sum(SPD'.*Zcmf)*Norm;
  [LD1(m,Lness) AD1(m,Lness) BD1(m,Lness)] = XYZtoLAB(X,Y,Z,Xw,Yw,Zw);  % convert to LAB
 end
end

nmin2 = floor((tmax2-1)/60);
LD2 = zeros(nmin2,LN,'double');
AD2 = zeros(nmin2,LN,'double');
BD2 = zeros(nmin2,LN,'double');

for Lness = 1:LN
 for m = 1:nmin2
  t = 60*(m-1)+1;                      % index to seconds
  r = uint8(255*RS2(t,Lness));         % get R,G,B digital values
  g = uint8(255*GS2(t,Lness));
  b = uint8(255*BS2(t,Lness));
  rgb = double([r g b]);
  SPD = zeros(s1count,1,'double');
  for k = 1:3
    SPD = SPD+squeeze(RGB_SPD(:,k))*RGBlut(k,rgb(k));  % make composite spectrum
  end
  X = sum(SPD'.*Xcmf)*Norm;
  Y = sum(SPD'.*Ycmf)*Norm;           % calculate tristimulus values
  Z = sum(SPD'.*Zcmf)*Norm;
  [LD2(m,Lness) AD2(m,Lness) BD2(m,Lness)] = XYZtoLAB(X,Y,Z,Xw,Yw,Zw);  % convert to LAB
 end
end

% Plot L*,a*,b* at one minute intervals

figure;  hold on;  grid on;  axis square;  rotate3d;
title('CIELAB a*-b* plane vs time for all lightness');
xlabel('a*');
ylabel('b*');
zlabel('L*');

for Lness = 1:LN
 h1 = plot3(AD1(:,Lness),BD1(:,Lness),LD1(:,Lness),':k','LineWidth',2);
 for m = 1:nmin1
  t = 60*(m-1)+1;                     % index to seconds
  r = RS1(t,Lness);                   % get R,G,B digital values
  g = GS1(t,Lness);
  b = BS1(t,Lness);
  rgb = double([r g b]);
  plot3(AD1(m,Lness),BD1(m,Lness),LD1(m,Lness),'ok','MarkerSize',10,...
    'MarkerEdgeColor','none','MarkerFaceColor',rgb);
 end
end

for Lness = 1:LN
 h2 = plot3(AD2(:,Lness),BD2(:,Lness),LD2(:,Lness),'-k','LineWidth',2);
 for m = 1:nmin2
  t = 60*(m-1)+1;                    % index to seconds
  r = RS2(t,Lness);                  % get R,G,B digital values
  g = GS2(t,Lness);
  b = BS2(t,Lness);
  rgb = double([r g b]);
  plot3(AD2(m,Lness),BD2(m,Lness),LD2(m,Lness),'ok','MarkerSize',10,...
      'MarkerEdgeColor','k','MarkerFaceColor',rgb);
 end
end

legend([h1 h2],name1,name2,'Location','SouthWest');

%% Compute LAB of match vs time from reconstructed spectra

Lness = 4;                           % selected lightness (index)

% Plot a*,b* as a function of time

figure;  hold on;  grid on;  axis square;  rotate3d;
title(sprintf('CIELAB a*-b* plane vs time, lightness = %d',90-5*Lness));
xlabel('a*');
ylabel('b*');
zlabel('Time (minutes)');
h1 = plot3(AD1(:,Lness),BD1(:,Lness),1:nmin1,':k','LineWidth',2);
for m = 1:nmin1
  t = 60*(m-1)+1;                     % index to seconds
  r = RS1(t,Lness);                   % get R,G,B digital values
  g = GS1(t,Lness);
  b = BS1(t,Lness);
  rgb = double([r g b]);
  plot3(AD1(m,Lness),BD1(m,Lness),m,'ok','MarkerSize',10,'MarkerEdgeColor','none',...
      'MarkerFaceColor',rgb);
end

h2 = plot3(AD2(:,Lness),BD2(:,Lness),1:nmin2,'-k','LineWidth',2);
for m = 1:nmin2
  t = 60*(m-1)+1;                    % index to seconds
  r = RS2(t,Lness);                  % get R,G,B digital values
  g = GS2(t,Lness);
  b = BS2(t,Lness);
  rgb = double([r g b]);
  plot3(AD2(m,Lness),BD2(m,Lness),m,'ok','MarkerSize',10,'MarkerEdgeColor','k',...
      'MarkerFaceColor',rgb);
end

legend([h1 h2],name1,name2,'Location','NorthWest');