%  Create a GUI panel with RGB colour sliders

function Colour_GUI

global red_slider;
global green_slider;
global blue_slider;

%  Create and hide the GUI as it is being constructed
  
f = figure('Visible','off','Position',[200,200,400,300]);
   
%  Construct the components
 
px = 50;  py = 100;  pw = 150;  ph = 150;
xsr = 250; xsg = 300; xsb = 350;
ysr = 100; ysg = 100; ysb = 100;
xsl =  20; ysl = 150;
tl = 30; tu = 10; xt = 24; yt = 16;
rgb = [0,0,0];
hred = uicontrol('Style','slider','Position',[xsr,ysr,xsl,ysl],...
                  'Callback',{@redslider_Callback});
hgreen = uicontrol('Style','slider','Position',[xsg,ysg,xsl,ysl],...
                  'Callback',{@greenslider_Callback});
hblue = uicontrol('Style','slider','Position',[xsb,ysb,xsl,ysl],...
                  'Callback',{@blueslider_Callback});
hr = uicontrol('Style','text','String','R','Position',[xsr,ysr-tl,xt,yt],...
                  'BackgroundColor',[0.8,0.5,0.5]);
hg = uicontrol('Style','text','String','G','Position',[xsg,ysr-tl,xt,yt],...
        'FontWeight','Bold','BackgroundColor',[0.5,0.8,0.5]);
hb = uicontrol('Style','text','String','B','Position',[xsb,ysr-tl,xt,yt],...
        'FontWeight','Bold','BackgroundColor',[0.5,0.5,0.8]);
red_slider = 0; green_slider = 0; blue_slider = 0;
UpdateColour();
align([hred,hgreen,hblue],'Top','None');
  
%  Make the GUI visible
  
set(f,'Visible','on');

% Helper functions

function redslider_Callback(hObject,eventdata,handles) 
  % Callback called when red colour slider is moved
  red_slider = get(hObject,'Value');
  UpdateColour();
end 

function greenslider_Callback(hObject,eventdata,handles) 
  % Callback called when green colour slider is moved
  green_slider = get(hObject,'Value');
  UpdateColour();
end

function blueslider_Callback(hObject,eventdata,handles) 
  % Callback called when blue colour slider is moved
  blue_slider = get(hObject,'Value');
  UpdateColour();
end

% Update RGB slider numeric values and colour panel

function UpdateColour()
  hrn = uicontrol('Style','text','String',sprintf('%3d',round(255*red_slider)),...
      'Position',[xsr-3,ysr+ysl+tu,xt,yt]);
  hrg = uicontrol('Style','text','String',sprintf('%3d',round(255*green_slider)),...
      'Position',[xsg-3,ysg+ysl+tu,xt,yt]);
  hrb = uicontrol('Style','text','String',sprintf('%3d',round(255*blue_slider)),...
      'Position',[xsb-3,ysb+ysl+tu,xt,yt]);
  rgb = [red_slider, green_slider, blue_slider];
  hp = uicontrol('Style','text','Position',[px,py,pw,ph],'BackgroundColor',rgb);
end

end

