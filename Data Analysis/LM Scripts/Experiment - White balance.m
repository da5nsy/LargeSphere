%----------------------------------------------------------------------
%  VISION SPHERE EXPERIMENT - WHITE BALANCE
%
%  Setup Phidget Interface via USB for following I/O:
%    Inputs from two sliders for A,B (analogue)
%    Input from push button (digital)
%
%----------------------------------------------------------------------

% Define data values

ON = 1;  OFF = 0;
TRUE = 1;  FALSE = 0;

pause on;                   % enable pausing

% Initialise Phidget controller

loadlibrary phidget21 phidget21Matlab.h;    % initialise Phidget library

ikptr = libpointer('int32Ptr',0);
er = calllib('phidget21','CPhidgetInterfaceKit_create',ikptr);  % set structure
if er ~= 0
    disp(strcat('Phidget allocation error_',sprintf('%d',er)));
end
ikhandle = get(ikptr, 'Value');

er = calllib('phidget21','CPhidget_open',ikhandle,-1);   % open device
if er ~= 0
    disp(strcat('Phidget open error_',sprintf('%d',er)));
end

er = calllib('phidget21','CPhidget_waitForAttachment',ikhandle,2500);
if er ~= 0
    disp(strcat('Phidget attachment error_',sprintf('%d',er)));
end

%% Set up for experiment

lambda = [400,420,440,460,480,500,520,540,560,580,600,620,640,660,680,700];
N = size(lambda,2);             % number of filters

% Define display panel for single white field

px = 800;  py = 400;  pw = 1000;  ph = 400;  % window location and dimensions
pw2 = floor(pw/2);  w = 200;
ax1 = pw2-w;  ax2 = pw2+w;                   % active area
f = figure('Visible','on','Position',[px,py,pw,ph]);
vbuf = zeros(ph,pw,3,'uint8');               % image buffer
vbuf(:,ax1:ax2,:) = 255;            % set white in central field
image(vbuf);                        % display target for setup

A_SLIDER = 1;                       % slider assignments for A,B
B_SLIDER = 0;
BUTTON = 1;                         % button input assignment
maxval = 1000;                      % maximum slider value
abuf = [0 0 0 0 0];                 % buffer for 5-point smoothing of slider
bbuf = [0 0 0 0 0];
paval = -1; pbval = -1;             % previous state of AB sliders
cfac = 40;                          % chroma scaling factor
lfac = 10;                          % lightness increment
rfac = 0.5;                         % scaling factor for random offset
LN = 5;                             % number of lightness levels per filter
dataptr = libpointer('int32Ptr',0);
LABmatch = zeros(3,LN,N,'double');  % slider values for match for each lamp
RGBmatch = zeros(3,LN,N,'double');  % converted R,G,B values for match for each lamp
Xw = 93.38; Yw = 100; Zw = 96.32;   % white tristimulus for D65

%% Main loop

for lam = 1:N                       % repeat for each wavelength of adapting light
 fprintf(1,'\n--> Filter wavelength = %d nm\n',lambda(lam));
 for n = 1:LN                       % repeat for each lightness level
  BUTTON_PRESSED = FALSE;           % reset state of push button
  azval = 0.5+rfac*(rand-0.5);      % randomise zero point on scale   
  bzval = 0.5+rfac*(rand-0.5);

% Get two slider settings and button status via Phidget controller

  while ~BUTTON_PRESSED
   er = calllib('phidget21','CPhidgetInterfaceKit_getSensorValue',ikhandle,A_SLIDER,dataptr);
   araw = get(dataptr,'Value');
   abuf(2:5) = abuf(1:4);           % shift buffer one place
   abuf(1) = double(araw);          % save new value
   aval = mean(abuf);
   er = calllib('phidget21','CPhidgetInterfaceKit_getSensorValue',ikhandle,B_SLIDER,dataptr);
   braw = get(dataptr,'Value');
   bbuf(2:5) = bbuf(1:4);           % shift buffer one place
   bbuf(1) = double(braw);          % save new value
   bval = mean(bbuf);
   er = calllib('phidget21','CPhidgetInterfaceKit_getInputState',ikhandle,BUTTON,dataptr);  
   butval = get(dataptr,'Value');
   if butval == ON
    BUTTON_PRESSED = TRUE;
   else
    L = 55+5*(LN+1-n); %lfac*(LN+1-n);                      % lightness
    as = double(maxval-aval)/maxval;        % slider value inverted (0 at top)
    bs = double(maxval-bval)/maxval;
    A = cfac*(azval-as)*(L/20);             % A relative to current centre
    B = cfac*(bzval-bs)*(L/20);             % B relative to current centre

% Convert slider settings to display R,G,B

    [X Y Z] = LABtoXYZ(L,A,B,Xw,Yw,Zw);  % convert LAB to XYZ
    [r g b gamutflag] = XYZtosRGBe(X/100,Y/100,Z/100); % convert to sRGB

% Update display field
 
    rgb = uint8(255*[r g b]);       % convert to integer 0-255
    vbuf(:,ax1:ax2,1) = rgb(1);     % RGB colour
    vbuf(:,ax1:ax2,2) = rgb(2);
    vbuf(:,ax1:ax2,3) = rgb(3);

    if (aval ~= paval) || (bval ~= pbval)           % has either slider moved?
      image(vbuf);                                  % yes, update image
%      fprintf(1,sprintf('L=%d,C=%d,H=%d RGB_raw = %4.2f,%4.2f,%4.2f  RGB_left = %3d,%3d,%3d  RGB_right = %3d,%3d,%3d\n',...
%        lval,cval,hval,r,g,b,rgbl,rgbr));          % and print out values
      paval = aval;  pbval = bval;                  % save current slider values
    end
    pause(0.1);                                     % allow time to complete
   end
  end

% Button pressed = best match achieved

  LABmatch(:,n,lam) = [L A B];           % save slider L,A,B values
  RGBmatch(:,n,lam) = [r g b];           % save converted R,G,B values
  paval = -1;                            % force display update
  pause(0.5);
 end

 fprintf(1,'\nMatching values for %d filter\n',lambda(lam));
 for n = 1:LN
  fprintf(1,'%d LAB= %5.2f,%6.2f,%5.2f  RGB_raw = %4.2f,%4.2f,%4.2f\n',...
      n,LABmatch(:,n,lam),RGBmatch(:,n,lam));
 end

% Flash display field to show ready for next filter

 zvbuf = vbuf;
 zvbuf(:,:,:) = 0;
 image(zvbuf);
 pause(0.5);
 image(vbuf);
 pause(0.5);
end

image(zvbuf);           % black display to show end
 
%% Plot results

for n = 1:LN
 figure;  hold on;
 title('LAB slider settings for white balance');
 L = LABmatch(1,n,1);
 text(420,L-5,sprintf('Lightness = %d',L));
 %axis([400 700 0 110]);
 ls = squeeze(LABmatch(1,n,:));
 plot(lambda,ls,'-k');                   % L
 hs = squeeze(LABmatch(2,n,:));
 plot(lambda,hs,'-r');                   % A
 cs = squeeze(LABmatch(3,n,:));
 plot(lambda,cs,'-b');                   % B
 xlabel('Filter wavelength (nm)');
 ylabel('Slider setting');
 legend('L','A','B','Location','SouthEast');
end

%% RGB

for n = 1:LN
 figure;  hold on;
 title('Normalised RGB values for white balance');
 %axis([400 700 -0.5 1.0]);
 L = LABmatch(1,n,1);
 rmax = max(RGBmatch(1,n,:));
 text(420,rmax,sprintf('Lightness = %d',L));
 rs = squeeze(RGBmatch(1,n,:));
 plot(lambda,rs,'-r');                   % R
 gs = squeeze(RGBmatch(2,n,:));
 plot(lambda,gs,'-g');                   % G
 bs = squeeze(RGBmatch(3,n,:));
 plot(lambda,bs,'-b');                   % B
 xlabel('Filter wavelength (nm)');
 ylabel('Normalised R,G,B');
 legend('R','G','B','Location','SouthEast');
end

%% Analyse differences across luminance levels

for i = 1:LN
  leg{:,i} = sprintf('L = %d',LABmatch(1,i,1));     % text for legend
end

% RGB

figure;  hold on;
title('Red display values relative to mean');
for n = 1:LN
  rm = geomean(RGBmatch(1,n,:));            % mean of R
  rs = squeeze(RGBmatch(1,n,:));            % data vector
  L = LABmatch(1,n,1)/100;
  cr = L+0.15;                              % colour scale blue-red
  cb = 1-cr;
  cg = 0.4;
  plot(lambda,rs/rm,'-','Color',[cr cg cb]);
end
xlabel('Filter wavelength (nm)');
ylabel('Normalised R');
%legend(leg,'Location','SouthEast');
plot([400 700],[1 1],':k');                 % dotted line at one

figure;  hold on;
title('Green display values relative to mean');
for n = 1:LN
  gm = geomean(RGBmatch(2,n,:));            % mean of G
  gs = squeeze(RGBmatch(2,n,:));            % data vector
  L = LABmatch(1,n,1)/100;
  cr = L+0.15;                              % colour scale blue-red
  cb = 1-cr;
  cg = 0.4;
  plot(lambda,gs/gm,'-','Color',[cr cg cb]);
end
xlabel('Filter wavelength (nm)');
ylabel('Normalised G');
%legend(leg,'Location','SouthWest');
plot([400 700],[1 1],':k');                 % dotted line at one

figure;  hold on;
title('Blue display values relative to mean');
for n = 1:LN
  bm = geomean(RGBmatch(3,n,:));            % mean of B
  bs = squeeze(RGBmatch(3,n,:));            % data vector
  L = LABmatch(1,n,1)/100;
  cr = L+0.15;                              % colour scale blue-red
  cb = 1-cr;
  cg = 0.4;
  plot(lambda,bs/bm,'-','Color',[cr cg cb]);
end
xlabel('Filter wavelength (nm)');
ylabel('Normalised B');
%legend(leg,'Location','SouthEast');
plot([400 700],[1 1],':k');             % dotted line at one

%% Plot results as 3D surfaces

figure;  hold on;  rotate3d;  grid on;
[XL,YL] = meshgrid(85:-5:10,lambda);    % lightness on X axis; wavelength on Y axis
ZL = squeeze(RGBmatch(1,:,:));          % red display values
for n = 1:LN
  rm = geomean(ZL(n,:));
  ZL(n,:) = ZL(n,:)/rm;                % normalise over mean of all wavelengths
end
surf(XL,YL,ZL');
xlabel('Lightness');
ylabel('Wavelength (nm)');
zlabel('Red display signal relative to mean');

figure;  hold on;  rotate3d;  grid on;
ZL = squeeze(RGBmatch(2,:,:));          % green display values
for n = 1:LN
  gm = geomean(ZL(n,:));
  ZL(n,:) = ZL(n,:)/gm;                % normalise over mean of all wavelengths
end
surf(XL,YL,ZL');
xlabel('Lightness');
ylabel('Wavelength (nm)');
zlabel('Green display signal relative to mean');

figure;  hold on;  rotate3d;  grid on;
ZL = squeeze(RGBmatch(3,:,:));          % blue display values
for n = 1:LN
  bm = geomean(ZL(n,:));
  ZL(n,:) = ZL(n,:)/bm;                % normalise over mean of all wavelengths
end
surf(XL,YL,ZL');
xlabel('Lightness');
ylabel('Wavelength (nm)');
zlabel('Blue display signal relative to mean');

%% LAB

figure;  hold on;
title('A slider values relative to mean');
for n = 1:LN
  am = mean(LABmatch(2,n,:));           % mean A
  as = squeeze(LABmatch(2,n,:));        % data vector
  c = 0.2*(LN-n);                       % grey value
  plot(lambda,as-am,'-','Color',[c c c]);
end
xlabel('Filter wavelength (nm)');
ylabel('Difference from mean A');
legend(leg,'Location','SouthEast');
plot([400 700],[0 0],':k');             % dotted line at zero

figure;  hold on;
title('B slider values relative to mean');
for n = 1:LN
  bm = mean(LABmatch(3,n,:));           % mean A
  bs = squeeze(LABmatch(3,n,:));        % data vector
  c = 0.2*(LN-n);                       % grey value
  plot(lambda,bs-bm,'-','Color',[c c c]);
end
xlabel('Filter wavelength (nm)');
ylabel('Difference from mean B');
legend(leg,'Location','NorthWest');
plot([400 700],[0 0],':k');             % dotted line at zero

%% Save data to text file

filename = fullfile('C:','Test','Experiment','White balance III.txt');
fp = fopen(filename,'w');
fprintf(fp,'Colour values for visually neutral field in sphere\n');

for lam = 1:N
  fprintf(fp,'\nFilter wavelength = %d nm\n',lambda(lam));
  for n = 1:LN
    fprintf(fp,'%d LAB = %5.2f,%6.2f,%5.2f  RGB = %6.3f,%6.3f,%6.3f\n',...
      n,LABmatch(:,n,lam),RGBmatch(:,n,lam));
 end
end

fclose(fp);
fprintf('Wrote data to file %s\n',filename);

save('White balance values III','LABmatch','RGBmatch');

%% Clean up

vbuf(:,:,:) = 0;
image(vbuf);                 % clear display
close(f);                    % close window

er = calllib('phidget21','CPhidget_close',ikhandle);  % release Phidget structure
if er ~= 0
  disp(strcat('Phidget close error ',sprintf('%d',er)));
end

er = calllib('phidget21','CPhidget_delete',ikhandle);
if er ~= 0
  disp(strcat('Phidget delete error ',sprintf('%d',er)));
end

disp('Bye bye');

%% Combine three datasets

load 'White balance values I'       % get first set of data (10,20,30,40,50)
RGB1 = RGBmatch; LAB1 = LABmatch;

load 'White balance values II'      % get second set of data (15,25,35,45,55)
RGB2 = RGBmatch; LAB2 = LABmatch;

load 'White balance values III'      % get second set of data (60,65,70,75,80,85)
RGB3 = RGBmatch; LAB3 = LABmatch;

LN = 16;
LABmatch = zeros(3,LN,N,'double');  % slider values for match for each lamp
RGBmatch = zeros(3,LN,N,'double');  % converted R,G,B values for match for each lamp

for n = 1:5
  RGBmatch(:,6+2*n,:) = RGB1(:,n,:);   % interleave two datasets
  RGBmatch(:,6+2*n-1,:) = RGB2(:,n,:);  
  LABmatch(:,6+2*n,:) = LAB1(:,n,:);   % interleave two datasets
  LABmatch(:,6+2*n-1,:) = LAB2(:,n,:);
end

for n = 1:6
  RGBmatch(:,n,:) = RGB3(:,n,:);       % append third dataset
  LABmatch(:,n,:) = LAB3(:,n,:);
end

save('White balance values combined','LABmatch','RGBmatch');
