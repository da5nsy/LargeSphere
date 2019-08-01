% Version of AveragingOverLStars that doesn't interpolate data.

% Let's try and edit this badboy to get some averages out.
% I reasoned in my ICVS presentation that we were safe with L* = 35 and 60
% (Lval = 6:11)


%% Plot 3D surfaces of LMSR vs wavelength and time for constant lightness

% lness = 6:11;              % lightness (fixed)
% WI = 400:20:700;
% wl = length(WI);
% LMSI = zeros(4,wl,TNM,'double');
% 
% for t = 1:TNM
%   for k = 1:4
%     v = mean(squeeze(LMSR(k,lness,t,:)));    % extract LMSR values for all wavelengths
%     LMSI(k,:,t) = v;  % interpolate to 5nm intervals
%   end
% end

LMSI=permute(squeeze(mean(LMSR(:,6:11,:,:),2)),[1,3,2]);

%% Plot figures
close all

%[XL,YL] = meshgrid(WI,1:TNM);

%figure;  hold on;  %rotate3d;  grid on;
%title('L cone excitation for wavelength vs time');
%ZL = squeeze(LMSI(1,:,:))';
%surf(XL,YL,ZL);
im = permute(LMSI(1,:,:),[3,2,1]);

%fakeim = cat(3,ones(30)/4,ones(30)/10,ones(30)/12);
%imagesc(fakeim)

imagesc(im)
set(gca,'YDir','normal')
%colorbar

xticks(1:16)
xticklabels(WI)

xlabel('Wavelength of adapting field (nm)');
ylabel('Time (min)');
%zlabel('L cone excitation');
%axis([WI(1) WI(wl) 1 TNM]);
    
% figure;  hold on;  rotate3d;  grid on;
% title('M cone excitation for wavelength vs time');
% ZL = squeeze(LMSI(2,:,:))';
% surf(XL,YL,ZL);
% xlabel('Wavelength of adapting field (nm)');
% ylabel('Time (min)');
% zlabel('M cone excitation');
% axis([WI(1) WI(wl) 1 TNM]);
% 
% figure;  hold on;  rotate3d;  grid on;
% title('S cone excitation for wavelength vs time');
% ZL = squeeze(LMSI(3,:,:))';
% surf(XL,YL,ZL);
% xlabel('Wavelength of adapting field (nm)');
% ylabel('Time (min)');
% zlabel('S cone excitation');
% axis([WI(1) WI(wl) 1 TNM]);
