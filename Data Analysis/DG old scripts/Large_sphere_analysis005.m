function Large_sphere_analysis005(obs,xbar2,ybar2,zbar2,lambdaCie2)
%005 corrects the Z tristimulus values rather than the xy chromaticities of
%the averages, multiplicitively, hopefully allowing a more full analysis to
%take place


%clc, clear, %close all

%obs='TR';   %'LM', 'TR', or 'DG'

%% Load CIE data
ciefile = fullfile('C:','Users','ucesars','Dropbox','UCL','Data',...
    'Colour Standards','CIE colorimetric data','CIE_colorimetric_tables.xls');

if exist('ciedata2','var')~=1
    ciedata2= xlsread(ciefile,'1931 col observer','A6:D86');
    lambdaCie2=ciedata2(:,1);
    xbar2=ciedata2(:,2);
    ybar2=ciedata2(:,3);
    zbar2=ciedata2(:,4);
    
    ciedata10= xlsread(ciefile,'1964 col observer','A6:D86');
    lambdaCie10=ciedata10(:,1);
    xbar10=ciedata10(:,2);
    ybar10=ciedata10(:,3);
    zbar10=ciedata10(:,4);
end

xbar2_101=interp1(lambdaCie2,xbar2,380:4:780,'spline');
ybar2_101=interp1(lambdaCie2,ybar2,380:4:780,'spline');
zbar2_101=interp1(lambdaCie2,zbar2,380:4:780,'spline');
%% Loads data

if strcmp(obs,'LM')
    rootdir = fullfile('C:','Users','ucesars','Dropbox','UCL','Data','Large Sphere','Lindsay time series - Sep 2012');
elseif strcmp(obs,'TR')
    rootdir = fullfile('C:','Users','ucesars','Dropbox','UCL','Data','Large Sphere','Tania  time series - Apr 2013');
elseif strcmp(obs,'DG')
    rootdir = fullfile('C:','Users','ucesars','Dropbox','UCL','Data','Large Sphere','Results - Oct 2016');
else display('No/incorrect observer entered')
end

cd(rootdir)
files= dir('*nm*.mat');

N = 10;                             % number of repetitions over time
LN = 16;                            % number of lightness levels per repeat

for j=1:length(files)
    
    load(fullfile(rootdir,files(j).name));  % load experimental results
    files(j).dataLAB=LABmatch;
    files(j).dataRGB=RGBmatch;
    files(j).Tmatch=Tmatch;
    
end

clear LABmatch RGBmatch RGBstart Tmatch j
%% Calibrate data

%load calibration file

if strcmp(obs,'LM')
    %%Same file as for TR, originally used as LM didn't remember doing one
    %%for his sessions, only for TR.
    %calFileLocation = fullfile('C:','Users','ucesars','Dropbox','UCL','Data','Large Sphere','Large LCD display measurement.mat');
    %%Later this file was found, which was created at the same time as the
    %%LM run, and so we assume it is appropriate to use.
    calFileLocation = fullfile('C:','Users','ucesars','Dropbox','UCL','Data','Large Sphere','Lindsay time series - Sep 2012','LCD display measurement.mat');
elseif strcmp(obs,'TR')
    calFileLocation = fullfile('C:','Users','ucesars','Dropbox','UCL','Data','Large Sphere','Large LCD display measurement.mat');
elseif strcmp(obs,'DG')
    calFileLocation = fullfile('C:','Users','ucesars','Dropbox','UCL','Data','Large Sphere','Large LCD display measurement - Oct 2016.mat');
else display('No/incorrect observer entered')
end
load(calFileLocation)

for trial=1:length(files)
    
    %For calibration based on XYZ vals
    %interpolate recorded values (sval) to required vals (0:1:255)
    XYZinterp=zeros(3,256,4);
    for i=1:3
        for j=1:4
            XYZinterp(i,:,j)=interp1(sval,XYZ(i,:,j),0:255,'spline');
        end
    end
    
    if strcmp(obs,'DG')
        needCorrecting=[440,460,500,540,560,580,600,620,660,680,700];
        if any(str2num(files(trial).name(1:3))==needCorrecting)
            XYZinterp(3,:,:)=XYZinterp(3,:,:)/1.7;
        end
        
    end
    
    %     %     % For calibration based on SPC
    %     %     SPCinterp=zeros(101,256,4);
    %     %     for i=1:101
    %     %         for j=1:4
    %     %             SPCinterp(i,:,j)=interp1(sval,Measurement(i,:,j),0:255,'spline');
    %     %         end
    %     %     end
    
    
    %Interp screen spectral data to same interval as CIE data
    RGBw_SPD = interp1(lambda,Measurement(:,21,4),lambdaCie2,'spline');
    RGBb_SPD = interp1(lambda,Measurement(:,1,4),lambdaCie2,'spline');
    
    Norm = 100/sum(RGBw_SPD.*ybar2);              % normalising factor
    
    DB = squeeze(RGBb_SPD);
    Xb = sum(RGBb_SPD.*xbar2)*Norm;
    Yb = sum(RGBb_SPD.*ybar2)*Norm;               % calculate white reference
    Zb = sum(RGBb_SPD.*zbar2)*Norm;
    %fprintf('%s Display black XYZ = %5.3f,%5.3f,%5.3f\n',files(trial).name(end-8:end-4),Xb,Yb,Zb);
    
    Xw = sum(RGBw_SPD.*xbar2)*Norm;
    Yw = sum(RGBw_SPD.*ybar2)*Norm;               % calculate white reference
    Zw = sum(RGBw_SPD.*zbar2)*Norm;
    %fprintf('%s Display white XYZ = %5.3f,%5.3f,%5.3f\n',files(trial).name(end-8:end-4),Xw,Yw,Zw);
    
    for j=1:3
        % Thresholding:
        % - original RGB values included out of gamut
        % selections, resulting in above 1 and below 0 values
        
        a=files(trial).dataRGB(j,:,:); %Create temporary variable
        a(files(trial).dataRGB(j,:,:)>1)=1; % Threshold to below 1
        a(files(trial).dataRGB(j,:,:)<0)=0; % Threshold to above 0
        files(trial).dataRGBcal(j,:,:)=uint8(a*255); %Rescale
        %dataRGBcal is not actually 'calibrated', perhaps innappropriate
        %naming
        
        files(trial).dataRGBcalgam(j,:,:)=files(trial).dataRGB(j,:,:)>1 ...
            |files(trial).dataRGB(j,:,:)<0;
        
    end
    
    files(trial).dataXYZ=zeros(3,LN,N);
    files(trial).dataxy=zeros(2,LN,N);
    files(trial).dataLABcal=zeros(3,LN,N);
    
    for j=1:LN
        for k=1:N
            %Using XYZinterp:
            files(trial).dataXYZ(:,j,k)=...
                (XYZinterp(:,files(trial).dataRGBcal(1,j,k)+1,1)...
                +XYZinterp(:,files(trial).dataRGBcal(2,j,k)+1,2)...
                +XYZinterp(:,files(trial).dataRGBcal(3,j,k)+1,3));
            
            %             R=files(trial).dataRGBcal(1,j,k);
            %             G=files(trial).dataRGBcal(2,j,k);
            %             B=files(trial).dataRGBcal(3,j,k);
            %             files(trial).combinedSpectrum(:,j,k)=SPCinterp(:,R+1,1)+SPCinterp(:,G+1,2)+SPCinterp(:,B+1,3);
            %             files(trial).dataXYZ(:,j,k)=[...
            %                 xbar2_101*files(trial).combinedSpectrum(:,j,k);
            %                 ybar2_101*files(trial).combinedSpectrum(:,j,k);
            %                 zbar2_101*files(trial).combinedSpectrum(:,j,k)];
            
            files(trial).dataxy(1,j,k)=...
                files(trial).dataXYZ(1,j,k)/sum(files(trial).dataXYZ(:,j,k));
            files(trial).dataxy(2,j,k)=...
                files(trial).dataXYZ(2,j,k)/sum(files(trial).dataXYZ(:,j,k));
            
            X=files(trial).dataXYZ(1,j,k);
            Y=files(trial).dataXYZ(2,j,k);
            Z=files(trial).dataXYZ(3,j,k);
            
            files(trial).datauv(1,j,k)=4*X/(X + 15*Y + 3*Z);
            files(trial).datauv(2,j,k)=9*Y/(X + 15*Y + 3*Z);
            
            files(trial).dataLABcal(:,j,k)=...
                XYZToLab(files(trial).dataXYZ(:,j,k),[Xw;Yw;Zw]);
        end
    end
end
%% Plot all data in LAB
figure('units','normalized','outerposition',[0 0 1 1])
hold on
axis('equal')
title(obs)
view(2)
LowL=30;    LL=18-LowL/5;
HighL=65;   HL=18-HighL/5;

xlim([-30 30])
ylim([-30 30])
zlim([LowL-5 HighL+5])

for i=1:length(files)
    %ignore un-corrected DG 480 data
    if strcmp(files(i).name,'480nm - time v2_3.mat')
    else
        
        LABdata=squeeze(mean(files(i).dataLABcal(:,HL:LL,:),3));
        s=plot3(LABdata(2,:),LABdata(3,:),LABdata(1,:));
        
        WL=(str2num(files(i).name(1:3)));
        LCol=xyz2rgb(ciedata2(ciedata2(:,1)==WL,2:4));
        LCol(LCol<0)=0;LCol(LCol>1)=1;
        s.Color=LCol;
        
        %drawnow
        %pause(.25)
        xlabel('a*')
        ylabel('b*')
        zlabel('L*')
        
        %LABdata=reshape(files(i).dataLABcal(:,HL:LL,:),3,(LL-HL+1)*N);
        %s=scatter3(LABdata(2,:),LABdata(3,:),LABdata(1,:),'r.');
        %s.CData=[(str2num(files(i).name(1:3))-400)/300,0,0];
        %s.SizeData=100;
        %text(LABdata(2),LABdata(3),files(i).name)
    end
end

% plot zero lines
currentaxes=gca;
plot([currentaxes.XLim],[0,0],'Color',[.8,.8,.8]);
plot([0,0],[currentaxes.YLim],'Color',[.8,.8,.8]);
%% Plot all data in u'v'

figure('DefaultAxesFontSize',24,'units','normalized','outerposition',[0 0 1 1])
hold on
axis('equal')
title(obs)

for trial=1:length(files)
    
    s(trial)=scatter(mean(reshape(files(trial).datauv(1,:,:),160,1)),...
              mean(reshape(files(trial).datauv(2,:,:),160,1)),...
              'k*');
    
%     WL=(str2num(files(trial).name(1:3)));
%     LCol=xyz2rgb(ciedata2(ciedata2(:,1)==WL,2:4));
%     LCol(LCol<0)=0;LCol(LCol>1)=1;
%     s.CData=LCol;
%     drawnow;
%     pause(0.1)
    
end
xlabel('u''')
ylabel('v''')
% PR650 Sphere Measurements

%arbitrary correction on XYZ data to get it to show up coloured on graph

%excluding high and low end spectra because of neutralish chromaticities
%due to low level = noise

% Load Sphere Spectra
load('Illumination in sphere');
spectralambda=[0 400:20:700];

% figure;  hold on;
% title('Spectral power distribution of filtered illumination in sphere');
% plot(lambda,spectra(:,1),'-y','LineWidth',2);
% plot(lambda,spectra(:,1),'-k');
% for i = 2:17
%   plot(lambda,spectra(:,i),'-k');
% end
% xlabel('Wavelength (nm)');
% ylabel('Radiant power');

% Calculate chromaticities and plot
for i=2:17 %5:14
    
    X(i)=xbar2_101*spectra(:,i)*1000;
    Y(i)=ybar2_101*spectra(:,i)*1000;
    Z(i)=zbar2_101*spectra(:,i)*1000;
    
    u_prime(i)= 4*X(i)/(X(i) + 15*Y(i) + 3*Z(i));
    v_prime(i)= 9*Y(i)/(X(i) + 15*Y(i) + 3*Z(i));
    
    s2=scatter(u_prime(i),v_prime(i),'filled','k');
%     CData=xyz2rgb([X(i),Y(i),Z(i)]);
%     CData(CData<0)=0;
%     s.CData=CData;
    text(u_prime(i)+0.01,v_prime(i),num2str(spectralambda(i)),'FontSize', 15)
end

for i=5:14%2:17
    plot([u_prime(i),s(i-1).XData] ,[v_prime(i),s(i-1).YData],'k');
end

axis('equal')

end
