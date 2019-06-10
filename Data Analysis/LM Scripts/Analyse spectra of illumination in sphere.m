%% Analyse spectra measured with PhotoResearch PR-650 spectroradiometer

load('Illumination in sphere');

% 'Spectral power distribution of light reflected from inner wall of sphere',...
% 'from a point just to the right of the circular aperture at rear.',...
% 'Light from Kodak S-AV 2000 projector with 150mm lens was filtered',...
% 'through 16 narrow-band filters and five broadband filters on lens.',...
% 'Measurements made with PR-650 spectroradiometer.',...
% 'Order: (no filter), 400, 420, 440, 460, 480, 500, 520, 540, 650, 580,',...
% '600, 620, 640, 660, 680, 700 nm, then red, orange, yellow, green, blue.',...
% 'Also spectrum of ambient room light, from overhead fluorescent luminaire',...
% 'reflected off white paper on table top.}';

%% Draw spectra of 16 narrow-band filters

figure;  hold on;
title('Spectral power distribution of filtered illumination in sphere');
plot(lambda,spectra(:,1),'-y','LineWidth',2);
plot(lambda,spectra(:,1),'-k');
for i = 2:17
  plot(lambda,spectra(:,i),'-k');
end
xlabel('Wavelength (nm)');
ylabel('Radiant power');

%% Draw spectra of broad-band filters

figure;  hold on;
title('Spectral power distribution of broadband filtered illumination');
plot(lambda,spectra(:,1),'-k');
plot(lambda,spectra(:,18),'-r');
plot(lambda,spectra(:,19),'-','Color',[1,0.7,0]);
plot(lambda,spectra(:,20),'-y');
plot(lambda,spectra(:,21),'-g');
plot(lambda,spectra(:,22),'-','Color',[0,0.7,1]);
xlabel('Wavelength (nm)');
ylabel('Radiant power');

%% Read standard observer CMFs - 380-800nm, 1nm intervals

cmffile = fullfile('C:','Research at UCL','Colour standards',...
            'CIE colorimetric data','StdObs-2deg-1nm.txt');
format = '%d %f %f %f';
fid = fopen(cmffile,'r');
[A,count] = fscanf(fid,format,[4,inf]);  % read the whole file into array A
fclose(fid);

len = length(lambda);           % force same length as TSR data
lambda1 = double(A(1,:));
lam = double(lambda);
Xcmf = interp1(lambda1,A(2,:),lam);
Ycmf = interp1(lambda1,A(3,:),lam);
Zcmf = interp1(lambda1,A(4,:),lam);

figure;
plot(lambda,Ycmf,'-k');

%% Calculate scaling factor for display luminance

load('LCD display measurement.mat');

swy = sum(RGBW_spectrum(:,4).*Ycmf(:));     % use Y as Vlambda
mfac = RGBW_XYZ(2,4)/swy;                   % white luminance

rlp = mfac*sum(RGBW_spectrum(:,1).*Ycmf);
glp = mfac*sum(RGBW_spectrum(:,2).*Ycmf);
blp = mfac*sum(RGBW_spectrum(:,3).*Ycmf);

rlm = RGBW_XYZ(2,1);
glm = RGBW_XYZ(2,2);
blm = RGBW_XYZ(2,3);

%% Calculate luminance of illuminated surround in sphere

spill = zeros(3,17,'double');       % XYZ of each spectrum

for n = 1:17
  spill(1,n) = mfac*sum(spectra(:,n).*Xcmf);
  spill(2,n) = mfac*sum(spectra(:,n).*Ycmf);
  spill(3,n) = mfac*sum(spectra(:,n).*Zcmf);
end

figure;  hold on;
title('Luminance of adapting field inside sphere');
plot(400:20:700,spill(2,2:17));
xlabel('Filter wavelength (nm)');
ylabel('Luminance (cd/m^2)');
