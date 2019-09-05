% Script to compute the CIELAB values for the adapting fields.

% We're really bending the rules of CIELab here:
% We're considering a single surface, lit by multiple different
% illuminants, and stating that the white point is a screen adjacent to
% this surface, which doesn't actually illuminate the surface at all.

% I don't feel great about this.

clear, clc, close all

% figure defaults
DGdisplaydefaults;

%% Load sphere measurements

load('C:\Users\cege-user\Dropbox\UCL\Data\LargeSphere\Hardware Data\Filter spectra\Illumination in sphere.mat','spectra','lambda')
lambda = double(lambda);
S_spectra = MakeItS(lambda); clear lambda
spectra = spectra(:,2:17);

%load('C:\Users\cege-user\Dropbox\UCL\Data\LargeSphere\Hardware Data\Calibration Data\LCD display measurement.mat','Measurement') %LM
load('C:\Users\cege-user\Dropbox\UCL\Data\LargeSphere\Hardware Data\Calibration Data\Large LCD display measurement.mat','Measurement') %TR
screen_spectra = Measurement(:,21,4);
clear Measurement
%% Load observer

load T_xyz1931.mat
T_xyz1931 = SplineCmf(S_xyz1931,T_xyz1931,S_spectra);

%% Calculate XYZ

XYZ = T_xyz1931 * spectra;
%Lab = XYZToLab(XYZ, [94.97,100,98.15]');
XYZ_screen = T_xyz1931 * screen_spectra;
XYZ_screenN = XYZ_screen/XYZ_screen(2)*100; %normalised
XYZN = XYZ/XYZ_screen(2)*100; %normalised, using the same factor as above
Lab = XYZToLab(XYZN, XYZ_screenN);

%% Plot

cols = jet(16);
% cols = XYZToSRGBPrimary(XYZ/max(XYZ(:)))';
% cols(cols>1) = 1;
% cols(cols<0) = 0;
adap = 400:20:700;

figure, hold on
plot3(Lab(2,:),Lab(3,:),Lab(1,:),'k','HandleVisibility','off')
for i = 1:16
    scatter3(Lab(2,i),Lab(3,i),Lab(1,i),[],cols(i,:),'filled','DisplayName',num2str(adap(i)))
end
grid off
xlabel('a*')
ylabel('b*')
zlabel('L*')
view(2)

legend('Location','bestoutside')

%save2pdf('C:\Users\cege-user\Dropbox\UCL\Data\LargeSphere\Data Analysis\figs\adapter1.pdf')

%%
figure, hold on
subplot(1,2,1), hold on
plot3(Lab(2,:),Lab(3,:),Lab(1,:),'k')
for i = 1:16
    scatter3(Lab(2,i),Lab(3,i),Lab(1,i),[],cols(i,:),'filled','DisplayName',num2str(adap(i)))
end
grid off
xlabel('a*')
ylabel('b*')
zlabel('L*')
view(0,0)

subplot(1,2,2), hold on
plot3(Lab(2,:),Lab(3,:),Lab(1,:),'k')
for i = 1:16
    scatter3(Lab(2,i),Lab(3,i),Lab(1,i),[],cols(i,:),'filled','DisplayName',num2str(adap(i)))
end
grid off
xlabel('a*')
ylabel('b*')
zticks([])
%axis equal
%cleanTicks
view(90,0)

%save2pdf('C:\Users\cege-user\Dropbox\UCL\Data\LargeSphere\Data Analysis\figs\adapter2.pdf')

%%

save('C:\Users\cege-user\Dropbox\UCL\Data\LargeSphere\Data Analysis\adaptingFieldCIELabTR','Lab')





