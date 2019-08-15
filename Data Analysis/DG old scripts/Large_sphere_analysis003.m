function Large_sphere_analysis003(obs)

clc, clear, %close all

obs='DG';   %'LM', 'TR', or 'DG'

%% Load CIE data
ciefile = fullfile('C:','Users','ucesars','Dropbox','UCL','Data',...
    'Colour Standards','CIE colorimetric data','CIE_colorimetric_tables.xls');

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
    

%     % For calibration based on SPC
%     SPCinterp=zeros(101,256,4);
%     for i=1:101
%         for j=1:4
%             SPCinterp(i,:,j)=interp1(sval,Measurement(i,:,j),0:255,'spline');
%         end
%     end
   
    
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
            
            files(trial).dataLABcal(:,j,k)=...
                XYZToLab(files(trial).dataXYZ(:,j,k),[Xw;Yw;Zw]);
        end
    end
end

% %% Plot all data in LAB
% 
% figure('units','normalized','outerposition',[0 0 1 1])
% hold on
% axis('equal')
% xlim([-60 40])
% ylim([-60 40])
% zlim([0 100])
% view(3)
% title(obs)
% 
% for i=1:length(files)
%     LABdata=reshape(files(i).dataLABcal,3,160);
%     s=scatter3(LABdata(2,:),LABdata(3,:),LABdata(1,:),'r.');
%     %s.CData=[LABdata(1,:)/100;1-(LABdata(1,:)/100);zeros(1,length(LABdata))]';
%     s.CData=[(str2num(files(i).name(1:3))-400)/300,0,0];
%     s.SizeData=100;
%     drawnow
%     pause(.25)
%     %s.CData=[0.2,0.2,0.2];
%     %s.SizeData=20;
%     xlabel('a*')
%     ylabel('b*')
%     zlabel('L*')
% end
% 
% % plot zero lines
% currentaxes=gca;
% plot([currentaxes.XLim],[0,0],'Color',[.8,.8,.8]);
% plot([0,0],[currentaxes.YLim],'Color',[.8,.8,.8]);
% 
% 
% % %% Plot all data in CIE xy
% % %plotCIE(3,xbar2,ybar2,zbar2)
% % hold on
% % axis('equal')
% % for i=1:length(files)
% %     xydata=reshape(files(i).dataxy,2,160);
% %     XYZdata=reshape(files(i).dataXYZ,3,160);
% %     scatter3(xydata(1,:),xydata(2,:),XYZdata(2,:)/100,'k.');
% % 
% % end
% % 
% % set(gca, 'DataAspectRatio', [repmat(min(diff(get(gca, 'XLim')), diff(get(gca, 'YLim'))), [1 2]) diff(get(gca, 'ZLim'))])


%% Plot all data in LAB

figure('units','normalized','outerposition',[0 0 1 1])
hold on
axis('equal')
xlim([-60 40])
ylim([-60 40])
title(obs)

for i=1:length(files)
    LABdata=mean(reshape(files(i).dataLABcal,3,160),2);
    s=scatter(LABdata(2,:),LABdata(3,:),'r.');
    s.CData=[(str2num(files(i).name(1:3))-400)/300,0,0];
    s.SizeData=100;
    drawnow
    %pause(.25)
    xlabel('a*')
    ylabel('b*')
    zlabel('L*')
    text(LABdata(2),LABdata(3),files(i).name)
end

% plot zero lines
currentaxes=gca;
plot([currentaxes.XLim],[0,0],'Color',[.8,.8,.8]);
plot([0,0],[currentaxes.YLim],'Color',[.8,.8,.8]);
