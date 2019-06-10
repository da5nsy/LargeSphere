%-------------------------------------------------------------------------
%  VISION SPHERE EXPERIMENT - WHITE BALANCE - TIME SERIES
%
%  Compare two sets of neutral match data
%
%-------------------------------------------------------------------------

w = 500;                            % wavelength of filter
N = 10;                             % number of repetitions over time
LN = 16;                            % number of lightness levels per repeat

dir1 = fullfile('C:','Research at UCL','Matlab scripts','Experiment',...
            'Results - Apr 2013');
fname1 = sprintf('%dnm - time v2.mat',w);   % load experimental results for wavelength
load(fullfile(dir1,fname1));
LAB1 = LABmatch;
RGB1 = RGBmatch;
T1 = Tmatch;

%dir2 = fullfile('C:','Research at UCL','Matlab scripts','Experiment',...
%        'Tania time series - Apr 2013');
fname2 = sprintf('%dnm - time.mat',w);   % load experimental results for wavelength
%load(fullfile(dir2,fname2));
load(fname2);
LAB2 = LABmatch;
RGB2 = RGBmatch;
T2 = Tmatch;

Lval = squeeze(LAB1(1,:,1,1));      % L values

%% Plot 3D surface of lightness vs iteration for RGB difference

RGBdiff = RGB1-RGB2;

rm = mean(mean(RGBdiff(1,3:LN-2,:)));  % mean ignoring two at each end of lightness range
gm = mean(mean(RGBdiff(2,3:LN-2,:)));
bm = mean(mean(RGBdiff(3,3:LN-2,:)));
fprintf('\nMean difference: R %6.4f  G %6.4f  B %6.4f\n',rm,gm,bm);
ram = mean(mean(abs(RGBdiff(1,3:LN-2,:))));
gam = mean(mean(abs(RGBdiff(2,3:LN-2,:))));
bam = mean(mean(abs(RGBdiff(3,3:LN-2,:))));
fprintf('\nMean absolute difference: R %6.4f  G %6.4f  B %6.4f\n',ram,gam,bam);

[XL,YL] = meshgrid(1:N,Lval);       % iteration on X axis; lightness on Y axis

figure;  hold on;  rotate3d;  grid on;
ZL = squeeze(RGBdiff(1,:,:));       % red difference values
surf(XL,YL,ZL);
title(sprintf('RGB difference, wavelength = %d',w));
xlabel('Iteration');
ylabel('Lightness');
zlabel('Red display value');
axis([1 N Lval(LN) Lval(1) min(min(ZL)) max(max(ZL))]);

figure;  hold on;  rotate3d;  grid on;
ZL = squeeze(RGBdiff(2,:,:));       % green difference values
surf(XL,YL,ZL);
title(sprintf('RGB difference, wavelength = %d',w));
xlabel('Iteration');
ylabel('Lightness');
zlabel('Green display value');
axis([1 N Lval(LN) Lval(1) min(min(ZL)) max(max(ZL))]);

figure;  hold on;  rotate3d;  grid on;
ZL = squeeze(RGBdiff(3,:,:));       % blue difference values
surf(XL,YL,ZL);
title(sprintf('RGB difference, wavelength = %d',w));
xlabel('Iteration');
ylabel('Lightness');
zlabel('Blue display value');
axis([1 N Lval(LN) Lval(1) min(min(ZL)) max(max(ZL))]);

%% Plot 3D surface of lightness vs iteration for relative RGB difference

RGBrel1 = zeros(3,LN,N,'double');
RGBrel2 = zeros(3,LN,N,'double');

for i = 1:LN
 for k = 1:3
   v = squeeze(RGB1(k,i,:));    % match values for one lightness, all iterations
   mv = mean(v);                % mean  over all iterations
   RGBrel1(k,i,:) = v/mv;       % normalise to average over all iterations
   v = squeeze(RGB2(k,i,:));    % match values for one lightness, all iterations
   mv = mean(v);                % mean  over all iterations
   RGBrel2(k,i,:) = v/mv;       % normalise to average over all iterations
 end
end

RGBdiff = RGBrel1-RGBrel2;

ram = 100*mean(mean(abs(RGBdiff(1,3:LN-2,:))));
gam = 100*mean(mean(abs(RGBdiff(2,3:LN-2,:))));
bam = 100*mean(mean(abs(RGBdiff(3,3:LN-2,:))));
fprintf('\nMean absolute relative difference (percent): R %4.2f  G %4.2f  B %4.2f\n',ram,gam,bam);

[XL,YL] = meshgrid(1:N,Lval(1:LN-2));       % iteration on X axis; lightness on Y axis

figure;  hold on;  rotate3d;  grid on;
ZL = squeeze(RGBdiff(1,1:LN-2,:));       % red difference values
surf(XL,YL,ZL);
title(sprintf('File #1: %s, file #2: %s, RGB difference, wavelength = %d',fname1,fname2,w));
xlabel('Iteration');
ylabel('Lightness');
zlabel('Red display value');
axis([1 N Lval(LN-2) Lval(1) min(min(ZL)) max(max(ZL))]);

figure;  hold on;  rotate3d;  grid on;
ZL = squeeze(RGBdiff(2,1:LN-2,:));       % green difference values
surf(XL,YL,ZL);
title(sprintf('File #1: %s, file #2: %s, RGB difference, wavelength = %d',fname1,fname2,w));
xlabel('Iteration');
ylabel('Lightness');
zlabel('Green display value');
axis([1 N Lval(LN-2) Lval(1) min(min(ZL)) max(max(ZL))]);

figure;  hold on;  rotate3d;  grid on;
ZL = squeeze(RGBdiff(3,1:LN-2,:));       % blue difference values
surf(XL,YL,ZL);
title(sprintf('File #1: %s, file #2: %s, RGB difference, wavelength = %d',fname1,fname2,w));
xlabel('Iteration');
ylabel('Lightness');
zlabel('Blue display value');
axis([1 N Lval(LN-2) Lval(1) min(min(ZL)) max(max(ZL))]);
   
%% Interpolate RGB to relative value every minute - first dataset

TM = T1;
t1h = squeeze(TM(1,1,1));               % time of first observation
t1m = squeeze(TM(2,1,1));
t1s = squeeze(TM(3,1,1));
t2h = squeeze(TM(1,LN,N));                 % time of last observation
t2m = squeeze(TM(2,LN,N));
t2s = squeeze(TM(3,LN,N));
tem1 = floor((3600*(t2h-t1h)+60*(t2m-t1m)+(t2s-t1s))/60);  % elapsed time for session (min)
RGBval1 = zeros(3,LN,tem1,'double');      % actual value per second (range 0-1)
RGBrel1 = zeros(3,LN,tem1,'double');      % actual value per second (range 0-1)
fprintf(1,'\nFile #1: %s, Start time %d:%d:%d  Session length = %d min\n',...
                    fname1,t1h,t1m,floor(t1s),tem1);

for k = 1:3
 for i = 1:LN
   v = squeeze(RGB1(k,i,:));      % match values for one lightness, all iterations
   th = squeeze(TM(1,i,:));       % hours
   tm = squeeze(TM(2,i,:));       % minutes
   ts = squeeze(TM(3,i,:));       % seconds
   te = floor((3600*(th-t1h)+60*(tm-t1m)+(ts-t1s))/60);  % elapsed time for each observation (sec)
   vi = interp1(te,v,1:tem1,'linear','extrap');  % interpolate to seconds  
   RGBval1(k,i,1:tem1) = vi;       % save interpolated lightness values
   vv = squeeze(RGB1(k,i,:));     % match values for one lightness, all iterations
   T = vv<0;  vv(T) = 0;          % clip negative values to zero
   mv = geomean(geomean(vv));     % mean  over all iterations
   RGBrel1(k,i,1:tem1) = vi/mv;    % normalise to average over all iterations
 end
end

%% Interpolate RGB to relative value every minute - second dataset

TM = T2;
t1h = squeeze(TM(1,1,1));               % time of first observation
t1m = squeeze(TM(2,1,1));
t1s = squeeze(TM(3,1,1));
t2h = squeeze(TM(1,LN,N));                 % time of last observation
t2m = squeeze(TM(2,LN,N));
t2s = squeeze(TM(3,LN,N));
tem2 = floor((3600*(t2h-t1h)+60*(t2m-t1m)+(t2s-t1s))/60);  % elapsed time for session (sec)
RGBval2 = zeros(3,LN,tem2,'double');      % actual value per second (range 0-1)
RGBrel2 = zeros(3,LN,tem2,'double');      % actual value per second (range 0-1)
fprintf(1,'\nFile #2: %s, Start time %d:%d:%d  Session length = %d min\n',...
                    fname2,t1h,t1m,floor(t1s),tem2);

for k = 1:3
 for i = 1:LN
   v = squeeze(RGB2(k,i,:));      % match values for one lightness, all iterations
   th = squeeze(TM(1,i,:));       % hours
   tm = squeeze(TM(2,i,:));       % minutes
   ts = squeeze(TM(3,i,:));       % seconds
   te = floor((3600*(th-t1h)+60*(tm-t1m)+(ts-t1s))/60);  % elapsed time for each observation (sec)
   vi = interp1(te,v,1:tem2,'linear','extrap');  % interpolate to seconds  
   RGBval2(k,i,1:tem2) = vi;       % save interpolated lightness values
   vv = squeeze(RGB2(k,i,:));     % match values for one lightness, all iterations
   T = vv<0;  vv(T) = 0;          % clip negative values to zero
   mv = geomean(geomean(vv));     % mean  over all iterations
   RGBrel2(k,i,1:tem2) = vi/mv;   % normalise to average over all iterations
 end
end

tem = min(tem1,tem2);              % length of shorter session

%% Plot 3D surface of lightness vs time for relative RGB difference

RGBrdiff = RGBrel1(:,:,1:tem)-RGBrel2(:,:,1:tem);

rm = mean(mean(RGBrdiff(1,3:LN-2,:)));  % mean ignoring two at each end of lightness range
gm = mean(mean(RGBrdiff(2,3:LN-2,:)));
bm = mean(mean(RGBrdiff(3,3:LN-2,:)));
fprintf('\nMean relative difference: R %6.4f  G %6.4f  B %6.4f\n',rm,gm,bm);
ram = mean(mean(abs(RGBrdiff(1,3:LN-2,:))));
gam = mean(mean(abs(RGBrdiff(2,3:LN-2,:))));
bam = mean(mean(abs(RGBrdiff(3,3:LN-2,:))));
fprintf('\nMean absolute difference: R %6.4f  G %6.4f  B %6.4f\n',ram,gam,bam);

[XL,YL] = meshgrid(1:tem,Lval);       % iteration on X axis; lightness on Y axis

figure;  hold on;  rotate3d;  grid on;
ZL = squeeze(RGBrdiff(1,:,:));       % red difference values
surf(XL,YL,ZL);
title(sprintf('RGB difference, wavelength = %d',w));
xlabel('Time (minutes)');
ylabel('Lightness');
zlabel('Red display value relative to lightness');
axis([1 tem Lval(LN) Lval(1) min(min(ZL)) max(max(ZL))]);

figure;  hold on;  rotate3d;  grid on;
ZL = squeeze(RGBrdiff(2,:,:));       % green difference values
surf(XL,YL,ZL);
title(sprintf('RGB difference, wavelength = %d',w));
xlabel('Time (minutes)');
ylabel('Lightness');
zlabel('Green display value relative to lightness');
axis([1 tem Lval(LN) Lval(1) min(min(ZL)) max(max(ZL))]);

figure;  hold on;  rotate3d;  grid on;
ZL = squeeze(RGBrdiff(3,:,:));       % blue difference values
surf(XL,YL,ZL);
title(sprintf('RGB difference, wavelength = %d',w));
xlabel('Time (minutes)');
ylabel('Lightness');
zlabel('Blue display value relative to lightness');
axis([1 tem Lval(LN) Lval(1) min(min(ZL)) max(max(ZL))]);

%% Interpolate LAB sliders to value every minute - first dataset

for k = 1:3
 for i = 1:LN
   v = squeeze(LAB1(k,i,:));      % match values for one lightness, all iterations
   th = squeeze(TM(1,i,:));       % hours
   tm = squeeze(TM(2,i,:));       % minutes
   ts = squeeze(TM(3,i,:));       % seconds
   te = floor((3600*(th-t1h)+60*(tm-t1m)+(ts-t1s))/60);  % elapsed time for each observation (sec)
   vi = interp1(te,v,1:tem1,'linear','extrap');  % interpolate to seconds  
   LABval1(k,i,1:tem1) = vi;       % save interpolated lightness values
 end
end

%% Interpolate LAB to relative value every minute - second dataset

for k = 1:3
 for i = 1:LN
   v = squeeze(LAB2(k,i,:));      % match values for one lightness, all iterations
   th = squeeze(TM(1,i,:));       % hours
   tm = squeeze(TM(2,i,:));       % minutes
   ts = squeeze(TM(3,i,:));       % seconds
   te = floor((3600*(th-t1h)+60*(tm-t1m)+(ts-t1s))/60);  % elapsed time for each observation (sec)
   vi = interp1(te,v,1:tem2,'linear','extrap');  % interpolate to seconds  
   LABval2(k,i,1:tem2) = vi;       % save interpolated lightness values
 end
end

%% Plot 3D surface of lightness vs time for relative LAB difference

LABdiff = LABval1(:,:,1:tem)-LABval2(:,:,1:tem);

am = mean(mean(LABdiff(2,:,:)));
bm = mean(mean(LABdiff(3,:,:)));
fprintf('\nMean difference: A %6.4f  B %6.4f\n',am,bm);

aam = mean(mean(abs(LABdiff(2,:,:))));
bam = mean(mean(abs(LABdiff(3,:,:))));
fprintf('\nMean absolute difference: A %6.4f  B %6.4f\n',aam,bam);

figure;  hold on;  rotate3d;  grid on;
ZL = squeeze(LABdiff(2,:,:));       % A difference values
surf(XL,YL,ZL);
title(sprintf('LAB difference, wavelength = %d',w));
xlabel('Time (minutes)');
ylabel('Lightness');
zlabel('A slider value');
axis([1 tem Lval(LN) Lval(1) min(min(ZL)) max(max(ZL))]);

figure;  hold on;  rotate3d;  grid on;
ZL = squeeze(LABdiff(3,:,:));       % B difference values
surf(XL,YL,ZL);
title(sprintf('LAB difference, wavelength = %d',w));
xlabel('Time (minutes)');
ylabel('Lightness');
zlabel('B slider value');
axis([1 tem Lval(LN) Lval(1) min(min(ZL)) max(max(ZL))]);
