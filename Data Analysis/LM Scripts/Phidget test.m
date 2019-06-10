% Test Phidget Interface Kit via USB
% The program toggles the two signals on digital outputs 0,1,2
% in order to switch LED test lamps on and off repeatedly.

% Define data values

ON = 1;
OFF = 0;
pause on;                   % enable pausing

% Initialise Phidget

loadlibrary phidget21 phidget21Matlab.h;    % initialise Phidget library

ikptr = libpointer('int32Ptr',0);
e = calllib('phidget21', 'CPhidgetInterfaceKit_create', ikptr);  % set structure
if e ~= 0
    disp(strcat('Phidget allocation error ',sprintf('%d',e)));
end
ikhandle = get(ikptr, 'Value');

e = calllib('phidget21', 'CPhidget_open', ikhandle, -1);   % open device
if e ~= 0
    disp(strcat('Phidget open error ',sprintf('%d',e)));
end
e = calllib('phidget21', 'CPhidget_waitForAttachment', ikhandle, 2500);
if e == 0
  for n = 1:100
      
    e = calllib('phidget21', 'CPhidgetInterfaceKit_setOutputState', ikhandle, 0, ON);
    if e ~= 0
      disp(strcat('Phidget ON error ',sprintf('%d',e)));
    end
    pause(0.5);
    e = calllib('phidget21', 'CPhidgetInterfaceKit_setOutputState', ikhandle, 0, OFF);
    if e ~= 0
      disp(strcat('Phidget OFF error ',sprintf('%d',e)));
    end
 
    e = calllib('phidget21', 'CPhidgetInterfaceKit_setOutputState', ikhandle, 1, ON);
    if e ~= 0
      disp(strcat('Phidget ON error ',sprintf('%d',e)));
    end
    pause(0.5);
    e = calllib('phidget21', 'CPhidgetInterfaceKit_setOutputState', ikhandle, 1, OFF);
    if e ~= 0
      disp(strcat('Phidget OFF error ',sprintf('%d',e)));
    end

    e = calllib('phidget21', 'CPhidgetInterfaceKit_setOutputState', ikhandle, 2, ON);
    if e ~= 0
      disp(strcat('Phidget ON error ',sprintf('%d',e)));
    end
    pause(0.5);
    e = calllib('phidget21', 'CPhidgetInterfaceKit_setOutputState', ikhandle, 2, OFF);
    if e ~= 0
      disp(strcat('Phidget OFF error ',sprintf('%d',e)));
    end

  end
else
  disp(strcat('Phidget attach error ',sprintf('%d',e)));
end

% Clean up

e = calllib('phidget21', 'CPhidget_close', ikhandle);  % release Phidget structure
if e ~= 0
  disp(strcat('Phidget close error ',sprintf('%d',e)));
end
 
e = calllib('phidget21', 'CPhidget_delete', ikhandle);
if e ~= 0
  disp(strcat('Phidget delete error ',sprintf('%d',e)));
end
 