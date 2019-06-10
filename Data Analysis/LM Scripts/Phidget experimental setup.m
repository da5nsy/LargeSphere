%----------------------------------------------------------------------
%
% Setup Phidget Interface via USB for following I/O:
%    Inputs from three sliders for R,G,B (analogue)
%    Input from push button (digital)
%    Output to eight LED lamps
%
%----------------------------------------------------------------------

% Define data values

ON = 1;  OFF = 0;
TRUE = 1;  FALSE = 0;

pause on;                   % enable pausing

% Initialise Phidget

loadlibrary phidget21 phidget21Matlab.h;    % initialise Phidget library

ikptr = libpointer('int32Ptr',0);
e = calllib('phidget21','CPhidgetInterfaceKit_create',ikptr);  % set structure
if e ~= 0
    disp(strcat('Phidget allocation error_',sprintf('%d',e)));
end
ikhandle = get(ikptr, 'Value');

e = calllib('phidget21','CPhidget_open',ikhandle,-1);   % open device
if e ~= 0
    disp(strcat('Phidget open error_',sprintf('%d',e)));
end

e = calllib('phidget21','CPhidget_waitForAttachment',ikhandle,2500);
if e ~= 0
    disp(strcat('Phidget attachment error_',sprintf('%d',e)));
end

%% Set up for experiment

lampcolour = {'White ','Red   ','Orange','Yellow',...
              'GreenY','GreenB','Blue  ','Violet'};

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

RED = 7;                    % slider assignments
GREEN = 6;
BLUE = 5;
BUTTON = 0;                 % button input assignment
maxval = 1000;              % maximum slider value
prval = -1; pgval = -1; pbval = -1;   % previous state of RGB sliders
zval = 0.15;                % zero point on scale
s = 1.5;                    % scaling factor for negative light
dataptr = libpointer('int32Ptr',0);
RGBmatch = zeros(3,8,'double'); % slider values for match for each lamp

%% Main loop

for n = 1:8
 BUTTON_PRESSED = FALSE;       % state of push button
 lampnum = n-1;
 e = calllib('phidget21','CPhidgetInterfaceKit_setOutputState',...
     ikhandle,lampnum,ON);

 while ~BUTTON_PRESSED
  e = calllib('phidget21','CPhidgetInterfaceKit_getSensorValue',ikhandle,RED,dataptr);
  rval = get(dataptr,'Value');
  e = calllib('phidget21','CPhidgetInterfaceKit_getSensorValue',ikhandle,GREEN,dataptr);
  gval = get(dataptr,'Value');
  e = calllib('phidget21','CPhidgetInterfaceKit_getSensorValue',ikhandle,BLUE,dataptr);
  bval = get(dataptr,'Value');
  e = calllib('phidget21','CPhidgetInterfaceKit_getInputState',ikhandle,BUTTON,dataptr);  
  butval = get(dataptr,'Value');
  if butval == ON
    BUTTON_PRESSED = TRUE;
  else
%    rs = double(maxval-rval)/maxval;    % sliders are inverted (0 at top)
%    gs = double(maxval-gval)/maxval;
%    bs = double(maxval-bval)/maxval;

    rs = double(rval)/maxval;            % scale raw slider value 0-1
    gs = double(gval)/maxval;
    bs = double(bval)/maxval;
    
    if (rs >= zval)
      rl = (rs-zval)/(1-zval);         % left side
      rr = 0;
    else
      rl = 0;                          % right side
      rr = s*(zval-rs);
    end
    
    if (gs >= zval)
      gl = (gs-zval)/(1-zval);         % left side
      gr = 0;
    else
      gl = 0;                          % right side
      gr = s*(zval-gs);
    end
    
    if (bs >= zval)
      bl = (bs-zval)/(1-zval);         % left side
      br = 0;
    else
      bl = 0;                          % right side
      br = s*(zval-bs);
    end
    
    rgbl = double([rl gl bl]);
    rgbr = double([rr gr br]);
    vbuf(:,ax1:ax2,1) = uint8(255*rl);     % RGB colour in left side
    vbuf(:,ax1:ax2,2) = uint8(255*gl);
    vbuf(:,ax1:ax2,3) = uint8(255*bl);
    vbuf(:,bx1:bx2,1) = uint8(255*rr);     % RGB colour in right side
    vbuf(:,bx1:bx2,2) = uint8(255*gr);
    vbuf(:,bx1:bx2,3) = uint8(255*br);
 %   fprintf(1,sprintf('R=%d, G=%d, B=%d, RGB = %4.2f,%4.2f,%4.2f\n',rval,gval,bval,rgb));
 %   pause(0.5);
    if (rval ~= prval) || (gval ~= pgval) || (bval ~= pbval)    % has any slider moved?
%      hp = uicontrol(f,'Style','text','Position',[1,1,pw,ph],'BackgroundColor',rgb);
       image(vbuf);
%      set(f,'CurrentAxes',ha)
%      axes(ha,'Visible','off');
       prval = rval;  pgval = gval;  pbval = bval;    % save current slider values
    end
    pause(0.1);
  end
 end
 
% Button pressed = best match achieved

 RGBmatch(:,n) = [rs gs bs];
 e = calllib('phidget21','CPhidgetInterfaceKit_setOutputState',...
     ikhandle,lampnum,OFF);
 pause(0.5);
end

%% Clean up

vbuf(:,:,:) = 0;
image(vbuf);                 % clear display
close(f);                    % close window

e = calllib('phidget21','CPhidget_close',ikhandle);  % release Phidget structure
if e ~= 0
  disp(strcat('Phidget close error ',sprintf('%d',e)));
end

e = calllib('phidget21','CPhidget_delete',ikhandle);
if e ~= 0
  disp(strcat('Phidget delete error ',sprintf('%d',e)));
end

fprintf(1,'\nRGB matching values\n');
for n = 1:8
  fprintf(1,'%n %s R=%5.3f R=%5.3f R=%5.3f\n',...
      n,lampcolour{n},RGBmatch(:,n));
end

disp('Bye bye');

%% ==== Flash all the lights ====

for i = 1:10
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
