
clc, clear, close all

% not working

%%
load('C:\Users\cege-user\Dropbox\UCL\Data\LargeSphere\Hardware Data\Filter spectra\Illumination in sphere.mat','spectra')
figure, plot(SToWls([380,4,101]),spectra(:,2:17))
axis tight
xlabel('Wavelength (nm)')
ylabel('Radiant Power')
%save2pdf('C:\Users\cege-user\Dropbox\UCL\Ongoing Work\Thesis\figs\LargeSphere\LSillum')

S_spectra = WlsToS([380:4:780]');


%%
load T_CIE_Y2.mat T_CIE_Y2 S_CIE_Y2

%figure,
%plot(SToWls(S_CIE_Y2),T_CIE_Y2)

%%

for i = 1:size(spectra,2)
   a(:,i) = SplineSpd(S_spectra,spectra(:,i),S_CIE_Y2);
   b(i) = T_CIE_Y2*a(:,i);
end

hold on
plot(SToWls(S_CIE_Y2),a(:,2:17))

%%

figure, plot(b)

