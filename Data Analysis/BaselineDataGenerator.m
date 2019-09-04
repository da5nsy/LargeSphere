function BaselineDataGenerator(A_SLIDER_val,B_SLIDER_val,azval,bzval)

%-------------------------------------------------------------------------
%  VISION SPHERE EXPERIMENT - WHITE BALANCE - TIME SERIES
%
% Alternative version where I can generate baseline data (no observer)
%
%-------------------------------------------------------------------------

% Define data values

TRUE = 1;  FALSE = 0;
JOYSTICK = 0;

pause on;                   % enable pausing

%% Set up for experiment

%dfile = fullfile('C:','Research at UCL','Experiment','Large LCD display measurement.mat');
%dfile = fullfile('C:','Research at UCL','Experiment','Filter Spectra','Large LCD display measurement - Oct 2016.mat');
dfile = fullfile('C:\Users\cege-user\Dropbox\UCL\Data\LargeSphere\Hardware Data\Calibration Data','Large LCD display measurement - Oct 2016.mat');
load(dfile);                        % load display data
%Xw = 99.04; Yw = 100; Zw = 151.30;  % normalised white tristimulus for display (cool, measured)
XYZw = XYZ(:,21,4);                  % white tristimulus values
Xw = 100*XYZw(1)/XYZw(2);
Yw = 100.0;
Zw = 100*XYZw(3)/XYZw(2);            % normalised white tristimulus for display

% Set up measurement panel on display

px = 295;  py = 409;  pw = 1000;  ph = 400;  % window location and dimensions
ph2 = floor(ph/2);  s = 8;
pw2 = floor(pw/2);  w = 200;
ax1 = pw2-w;  ax2 = pw2+w;              % active area
fh = figure('Visible','on','Position',[px,py,pw,ph]);  % create figure
set(fh,'Color','k','MenuBar','none');   % set figure background to black
vbuf = zeros(ph,pw,3,'uint8');          % image buffer
dbuf = zeros(ph,ph,3,'uint8');          % display template
rbuf = zeros(ph,ph,3,'uint8');
tbuf = zeros(ph,ph,3,'uint8');          % circular template
rad = 100;                              % radius of circular target
for i = 1:pw
  for j = 1:ph
    r = sqrt((ph2-i)^2+(ph2-j)^2);      % distance from centre
    if (r<=rad) 
        tbuf(j,i,:) = 1; 
    end    % 1 inside circle
  end
end
Tcircle = (tbuf==1);                    % logical 'true' inside circle

dbuf(Tcircle) = 255;                    % set white in central field
dbuf(ph2-s:ph2+s,ph2-s:ph2+s,:) = 0;    % small black square at centre
vbuf(:,ax1:ax2-1,:) = dbuf;             % fill only circle
image(vbuf);                            % display target for setup

% Define control channels

if JOYSTICK
  A_SLIDER = 0;                     % joystick assignments for A,B
  B_SLIDER = 1;
else
  A_SLIDER = 1;                     % slider assignments for A,B
  B_SLIDER = 0;
end
BUTTON = 1;                         % button input assignment

% Define constants for interaction

ur =  25;  uy = 92;                 % CIELAB unique hues for R,Y,G,B (degrees)
ug = 163;  ub = 253;                % (Derefeldt et al, 2007, CR&A 11(2)148-152)
sur = sind(ur);  cur = cosd(ur);
suy = sind(uy);  cuy = cosd(uy);
sug = sind(ug);  cug = cosd(ug);
sub = sind(ub);  cub = cosd(ub);

maxval = 1000;                      % maximum slider value
abuf = [0 0 0 0 0];                 % buffer for 5-point smoothing of slider
bbuf = [0 0 0 0 0];
paval = -1; pbval = -1;             % initialise 'previous' state of AB sliders
cfac = 50;                          % chroma scaling factor
lfac = 5;                           % lightness increment
rfac = 0.5;                         % scaling factor for random offset

%% Main loop

N = 1;                             % number of repetitions over time
LN = 16;                            % number of lightness levels per repeat
LABmatch = zeros(3,LN,N,'double');  % slider values for match for each lamp
RGBmatch = zeros(3,LN,N,'double');  % converted R,G,B values for match for each lamp
Tmatch = zeros(3,LN,N,'double');    % time (H,M,S) for match

for t = 1:N                         % repeat in time
 fprintf(1,'\n--> Repetition count = %d\n',t);
 for n = 1:LN                       % repeat for each lightness level
  BUTTON_PRESSED = FALSE;           % reset state of push button
  butval = 0;
%   azval = 0.5+rfac*(rand-0.5);      % randomise zero point on scale   
%   bzval = 0.5+rfac*(rand-0.5);

% Get two slider settings and button status via Phidget controller

  while ~BUTTON_PRESSED   
   araw = A_SLIDER_val;
   abuf(2:5) = abuf(1:4);           % shift buffer one place
   abuf(1) = double(araw);          % save new value
   aval = mean(abuf);
   braw = B_SLIDER_val;
   bbuf(2:5) = bbuf(1:4);           % shift buffer one place
   bbuf(1) = double(braw);          % save new value
   bval = mean(bbuf);
   butval = butval + 1;
   if butval > 15
    BUTTON_PRESSED = TRUE;
   else
    L = 85-lfac*(n-1);              % lightness (downward sequence)
    as = double(aval)/maxval;       % slider value inverted (0 at top)
    bs = double(bval)/maxval;
    az = as-azval;                  % A relative to current centre
    bz = bs-bzval;                  % B relative to current centre
    if JOYSTICK
      az = -az;
    end

% Convert slider settings to display R,G,B

    if az>0
      au = az*cur;                 % movement towards red
      av = az*sur;
    else
      au = -az*cug;                % movement towards green
      av = -az*sug;
    end
    if bz>0
      bu = bz*cuy;                 % movement towards yellow
      bv = bz*suy;
    else
      bu = -bz*cub;                % movement towards blue
      bv = -bz*sub;
    end
    A = cfac*(au+bu);  B = cfac*(av+bv);   % sum components
    A = A*(L/20);  B = B*(L/20);          % adjust for L level
    [X, Y, Z] = LABtoXYZ(L,A,B,Xw,Yw,Zw);    % convert LAB to XYZ
    [r, g, b, gamutflag] = XYZtosRGBe(X/100,Y/100,Z/100); % convert to sRGB

% Update display field

    rr = uint8(255*r);
    gg = uint8(255*g);
    bb = uint8(255*b);
    rbuf(:,:,1) = rr;         % fill buffer with same colour
    rbuf(:,:,2) = gg;
    rbuf(:,:,3) = bb;
    dbuf(Tcircle) = rbuf(Tcircle);      % set colour in central field
    vbuf(:,ax1:ax2-1,:) = dbuf;         % update circular target
    
    image(vbuf);                                  % yes, update image
    
    paval = aval;  pbval = bval;                  % save current slider values
    %pause(0.05);                                    % allow time to complete
   end
  end

% Button pressed = best match achieved

  LABmatch(:,n,t) = [L A B];           % save slider L,A,B values
  RGBmatch(:,n,t) = [r g b];           % save converted R,G,B values
  cl = clock;                          % get current time
  Tmatch(:,n,t) = cl(4:6);             % save H,M,S
  paval = -1;                          % force display update
  pause(0.01);                            % wait a second to ensure button debounce
 end

 fprintf(1,'\nMatching values for %d filter\n');
 for n = 1:LN
  fprintf(1,'%d LAB= %5.2f,%6.2f,%5.2f  RGB_raw = %4.2f,%4.2f,%4.2f\n',...
      n,LABmatch(:,n,t),RGBmatch(:,n,t));
 end
end

% Finishing flourish on display

for n = 1:10
  vbuf(:,:,:) = floor(150*n/10);
  image(vbuf);           % ramp display up
  pause(0.01);
end
for n = 1:10
  vbuf(:,:,:) = floor(150*(10-n)/10);
  image(vbuf);           % ramp display down
  pause(0.01);
end

%% Save data to text and binary files
filename = ['C:\Users\cege-user\Dropbox\UCL\Data\LargeSphere\Experimental Data\BaselineData\',num2str(A_SLIDER_val),'-',num2str(B_SLIDER_val),'-',num2str(azval*100),'-',num2str(bzval*100)]; %*100 to avoid decimal places in filenames #hacky
fp = fopen([filename,'.txt'],'w');
fprintf(fp,'Colour values for visually neutral field in sphere\n');
fprintf(fp,'Date %d-%d-%d\n',cl(3),cl(2),cl(1));
fprintf(fp,'\nWavelength = %d nm\n');
%fprintf(fp,'\nNo filter (Black)');

for t = 1:N
  fprintf(fp,'\nRepetition = %d\n',t);
  for n = 1:LN
    tt = squeeze(Tmatch(:,n,t));
    th = tt(1);  tm = tt(2);  ts = round(tt(3));
    fprintf(fp,'%2d Time %2d:%02d:%02d  LAB = %5.2f,%6.2f,%5.2f  RGB = %6.3f,%6.3f,%6.3f\n',...
      n,th,tm,ts,LABmatch(:,n,t),RGBmatch(:,n,t));
 end
end

fclose(fp);
fprintf(1,'Wrote data to file %s\n',filename);

Wref = [Xw Yw Zw];
save(filename,'Tmatch','LABmatch','RGBmatch','Wref');  % save file

%% End session and close Phidget interface

disp('Bye bye');

end
