%----------------------------------------------------------------------
%  VISION SPHERE EXPERIMENT
%
%  Setup Phidget Interface via USB for following I/O:
%    Inputs from three sliders for L,H,C (analogue)
%    Input from push button (digital)
%    Output to eight LED lamps
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

lampcolour = {'white','red','orange','yellow',...
              'greenY','greenB','blue','magenta'};

lrange = [2 4 7];               % range of lamps in use
lambda = [400,420,440,460,480,500,520,540,560,580,600,620,640,660,680,700];
N = size(lambda,2);             % number of filters

% Define display panel for bipartite colour fields

px = 100;  py = 400;  pw = 1600;  ph = 400;  % window dimensions
pw2 = floor(pw/2);  w = 200;
ax1 = pw2-w;  ax2 = pw2-1;                   % active area (left)
bx1 = pw2+2;  bx2 = pw2+w;                   % active area (right)
f = figure('Visible','on','Position',[px,py,pw,ph]);
vbuf = zeros(ph,pw,3,'uint8');               % image buffer
vbuf(:,ax1:ax2,1) = 255;            % red in left side
vbuf(:,bx1:bx2,2) = 255;            % green in right side
image(vbuf);                        % display target for setup

L_SLIDER = 5;                     % slider assignments for L,C,H
C_SLIDER = 6;
H_SLIDER = 7;
BUTTON = 0;                 % button input assignment
maxval = 1000;              % maximum slider value
plval = -1; phval = -1; pcval = -1;   % previous state of RGB sliders
zval = 0.15;                % zero point on scale
s = 1;                    % scaling factor for negative R,G,B (out of gamut)
dataptr = libpointer('int32Ptr',0);
LCHmatch = zeros(3,8,N,'double');  % slider values for match for each lamp
RGBmatch = zeros(3,8,N,'double');  % converted R,G,B values for match for each lamp
Xw = 93.38; Yw = 100; Zw = 96.32;  % white tristimulus for D65
r0 = 0;  g0 = 0;  b0 = 0;    % grey offset behind LED stimulus (normalised RGB)

%% Main loop

for lam = 1:N                   % repeat for each wavelength of adapting light
 fprintf(1,'\n--> Filter wavelength = %d nm\n',lambda(lam));
 for n = lrange                 % repeat for each LED target
  BUTTON_PRESSED = FALSE;       % reset state of push button
  lampnum = n-1;
  er = calllib('phidget21','CPhidgetInterfaceKit_setOutputState',ikhandle,lampnum,ON);

  while ~BUTTON_PRESSED

% Get three slider settings and button status via Phidget controller

   er = calllib('phidget21','CPhidgetInterfaceKit_getSensorValue',ikhandle,L_SLIDER,dataptr);
   lval = get(dataptr,'Value');
   er = calllib('phidget21','CPhidgetInterfaceKit_getSensorValue',ikhandle,C_SLIDER,dataptr);
   cval = get(dataptr,'Value');
   er = calllib('phidget21','CPhidgetInterfaceKit_getSensorValue',ikhandle,H_SLIDER,dataptr);
   hval = get(dataptr,'Value');
   er = calllib('phidget21','CPhidgetInterfaceKit_getInputState',ikhandle,BUTTON,dataptr);  
   butval = get(dataptr,'Value');
   if butval == ON
    BUTTON_PRESSED = TRUE;
   else
%   ls = double(maxval-lval)/maxval;     % sliders are inverted (0 at top)
%   cs = double(maxval-cval)/maxval;
%   hs = double(maxval-hval)/maxval;

    ls = 100*double(lval)/maxval;        % scale raw slider values
    cs = 120*double(cval)/maxval;
    hs = 360*double(hval)/maxval;
    as = cs*cosd(hs);                    % convert polar (CH) to cartesian (AB)
    bs = cs*sind(hs);

% Convert slider settings to display R,G,B

    [X Y Z] = LABtoXYZ(ls,as,bs,Xw,Yw,Zw);  % convert LAB to XYZ
    [r g b gamutflag] = XYZtosRGBe(X/100,Y/100,Z/100); % convert to sRGB

    if (r >= 0)
      rl = r;                       % left side
      rr = r0;
    else
      rl = 0;                       % right side
      rr = r0-s*r;
    end

    if (g >= 0)
      gl = g;                       % left side
      gr = g0;
    else
      gl = 0;                       % right side
      gr = g0-s*g;
    end

    if (b >= 0)
      bl = b;                       % left side
      br = b0;
    else
      bl = 0;                       % right side
      br = b0-s*b;
    end

% Update bipartite fields
 
    rgbl = uint8(255*[rl gl bl]);    % convert to integer 0-255
    rgbr = uint8(255*[rr gr br]);
    vbuf(:,ax1:ax2,1) = rgbl(1);     % RGB colour on left side
    vbuf(:,ax1:ax2,2) = rgbl(2);
    vbuf(:,ax1:ax2,3) = rgbl(3);
    vbuf(:,bx1:bx2,1) = rgbr(1);     % RGB colour on right side
    vbuf(:,bx1:bx2,2) = rgbr(2);
    vbuf(:,bx1:bx2,3) = rgbr(3);

    if (lval ~= plval) || (hval ~= phval) || (cval ~= pcval)    % has any slider moved?
      image(vbuf);                                  % yes, update image
%      fprintf(1,sprintf('L=%d,C=%d,H=%d RGB_raw = %4.2f,%4.2f,%4.2f  RGB_left = %3d,%3d,%3d  RGB_right = %3d,%3d,%3d\n',...
%        lval,cval,hval,r,g,b,rgbl,rgbr));           % and print out values
      plval = lval;  phval = hval;  pcval = cval;   % save current slider values
    end
    pause(0.1);                                     % allow time to complete
   end
  end

% Button pressed = best match achieved

  LCHmatch(:,n,lam) = [ls hs cs];       % save slider L,C,H values
  RGBmatch(:,n,lam) = [r g b];          % save converted R,G,B values
  er = calllib('phidget21','CPhidgetInterfaceKit_setOutputState',ikhandle,lampnum,OFF);
  pause(0.5);
 end

 fprintf(1,'\nMatching values for %d filter\n',lambda(lam));
 for n = lrange
  fprintf(1,'%d %s  LHC = %5.2f,%6.2f,%5.2f  RGB_raw = %4.2f,%4.2f,%4.2f\n',...
      n,lampcolour{n},LCHmatch(:,n,lam),RGBmatch(:,n,lam));
 end

% Flash lights to show ready for next filter

 for lampnum = 6:-1:1
  er = calllib('phidget21','CPhidgetInterfaceKit_setOutputState',ikhandle,lampnum,ON);
  pause(0.2);
  er = calllib('phidget21','CPhidgetInterfaceKit_setOutputState',ikhandle,lampnum,OFF);
  pause(0.2);
 end
end

%% Plot results

for n = lrange
 figure;  hold on;
 title(sprintf('LCH slider settings for matching %s LED',lampcolour{n}));
 %axis([400 700 0 110]);
 ls = squeeze(LCHmatch(1,n,:));
 plot(lambda,ls,'-k');                   % L
 hs = squeeze(LCHmatch(2,n,:));
 plot(lambda,hs,'-b');                   % H
 cs = squeeze(LCHmatch(3,n,:));
 plot(lambda,cs,'-r');                   % C
 xlabel('Filter wavelength (nm)');
 ylabel('Slider setting');
 legend('L','H','C','Location','East');
end

for n = lrange
 figure;  hold on;
 title(sprintf('Normalised RGB values for matching %s LED',lampcolour{n}));
 %axis([400 700 -0.5 1.0]);
 rs = squeeze(RGBmatch(1,n,:));
 plot(lambda,rs,'-r');                   % R
 gs = squeeze(RGBmatch(2,n,:));
 plot(lambda,gs,'-g');                   % G
 bs = squeeze(RGBmatch(3,n,:));
 plot(lambda,bs,'-b');                   % B
 xlabel('Filter wavelength (nm)');
 ylabel('Normalised R,G,B');
 legend('R','G','B','Location','East');
end

%% Analyse differences across lamps

col = [ 0,0,0; 1,0,0; 1,0.5,0; 1,0.9,0; 0,1,0; 0,1,0.7; 0,0,1; 1,0,1 ];

% RGB

figure;  hold on;
title('R slider values relative to mean');
plot([400 700],[0 0],':k');             % dotted line at zero
for n = lrange
  rm = mean(RGBmatch(1,n,:));           % mean L
  rs = squeeze(RGBmatch(1,n,:));        % data vector
  plot(lambda,rs-rm,'-','Color',col(n,:));
end
xlabel('Filter wavelength (nm)');
ylabel('Difference from mean R');

figure;  hold on;
title('G slider values relative to mean');
plot([400 700],[0 0],':k');             % dotted line at zero
for n = 1:8
  gm = mean(RGBmatch(2,n,:));           % mean G
  gs = squeeze(RGBmatch(2,n,:));        % data vector
  plot(lambda,gs-gm,'-','Color',col(n,:));
end
xlabel('Filter wavelength (nm)');
ylabel('Difference from mean G');

figure;  hold on;
title('B slider values relative to mean');
plot([400 700],[0 0],':k');             % dotted line at zero
for n = lrange
  bm = mean(RGBmatch(3,n,:));           % mean B
  bs = squeeze(RGBmatch(3,n,:));        % data vector
  plot(lambda,bs-bm,'-','Color',col(n,:));
end
xlabel('Filter wavelength (nm)');
ylabel('Difference from mean B');

%% LCH

figure;  hold on;
title('L slider values relative to mean');
plot([400 700],[0 0],':k');             % dotted line at zero
for n = lrange
  lm = mean(LCHmatch(1,n,:));           % mean L
  ls = squeeze(LCHmatch(1,n,:));        % data vector
  plot(lambda,ls-lm,'-','Color',col(n,:));
end
xlabel('Filter wavelength (nm)');
ylabel('Difference from mean L');

figure;  hold on;
title('H slider values relative to mean');
plot([400 700],[0 0],':k');             % dotted line at zero
for n = lrange
  hm = mean(LCHmatch(2,n,:));           % mean L
  hs = squeeze(LCHmatch(2,n,:));        % data vector
  plot(lambda,hs-hm,'-','Color',col(n,:));
end
xlabel('Filter wavelength (nm)');
ylabel('Difference from mean H');

figure;  hold on;
title('C slider values relative to mean');
plot([400 700],[0 0],':k');             % dotted line at zero
for n = lrange
  cm = mean(LCHmatch(3,n,:));           % mean L
  cs = squeeze(LCHmatch(3,n,:));        % data vector
  plot(lambda,cs-cm,'-','Color',col(n,:));
end
xlabel('Filter wavelength (nm)');
ylabel('Difference from mean C');

%% Save data to text file

filename = fullfile('C:','Test','Experiment','LM results - red,yellow,blue 10-Aug.txt');
fp = fopen(filename,'w');
fprintf(fp,'Colour matching values for bipartite field in sphere\n',lampcolour{n});

for lam = 1:N
  fprintf(fp,'\nFilter wavelength = %d nm\n',lambda(lam));
  for n = lrange
    fprintf(fp,'%d LHC = %5.2f,%6.2f,%5.2f  RGB_raw = %6.3f,%6.3f,%6.3f  %s\n',...
      n,LCHmatch(:,n,lam),RGBmatch(:,n,lam),lampcolour{n});
 end
end

fclose(fp);
fprintf('Wrote data to file %s\n',filename);

save('Match values','LCHmatch','RGBmatch');

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

%% ==== Flash all the lights ====

for i = 1:5
 for n = 0:7
  e = calllib('phidget21','CPhidgetInterfaceKit_setOutputState',ikhandle,n,ON);
  if e ~= 0
    disp(strcat('Phidget ON error ',sprintf('%d',e)));
  end
  pause(1);
  e = calllib('phidget21','CPhidgetInterfaceKit_setOutputState',ikhandle,n,OFF);
  if e ~= 0
    disp(strcat('Phidget OFF error ',sprintf('%d',e)));
  end
  pause(0.1);
 end
end
