%----------------------------------------------------------------------
%
% Setup Phidget Interface via USB for LUMINAR dual-wavelength lights
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
else
    disp('Phidget interface open');
end

e = calllib('phidget21','CPhidget_waitForAttachment',ikhandle,2500);
if e ~= 0
    disp(strcat('Phidget attachment error_',sprintf('%d',e)));
else
    disp('Phidget attached and ready');    
end

%% Flash lights ten times

BLUE = 0;                   % assignment for violet 405nm
RED = 1;                    % assignment for infrared 850nm

for n = 1:3
   e = calllib('phidget21','CPhidgetInterfaceKit_setOutputState',...
     ikhandle,BLUE,ON);
   if e ~= 0
    disp(strcat('BLUE ON error ',sprintf('%d',e)));
   end
   pause(0.5);
   e = calllib('phidget21','CPhidgetInterfaceKit_setOutputState',...
     ikhandle,BLUE,OFF);
   if e ~= 0
    disp(strcat('BLUE OFF error ',sprintf('%d',e)));
   end
   pause(1);
   e = calllib('phidget21','CPhidgetInterfaceKit_setOutputState',...
     ikhandle,RED,ON);
   if e ~= 0
    disp(strcat('RED ON error ',sprintf('%d',e)));
   end
   pause(0.5);
   e = calllib('phidget21','CPhidgetInterfaceKit_setOutputState',...
     ikhandle,RED,OFF);
   if e ~= 0
    disp(strcat('RED OFF error ',sprintf('%d',e)));
   end
   pause(1);
end

%% Clean up

e = calllib('phidget21','CPhidget_close',ikhandle);  % release Phidget structure
if e ~= 0
  disp(strcat('Phidget close error ',sprintf('%d',e)));
end

e = calllib('phidget21','CPhidget_delete',ikhandle);
if e ~= 0
  disp(strcat('Phidget delete error ',sprintf('%d',e)));
end

disp('Bye bye');
