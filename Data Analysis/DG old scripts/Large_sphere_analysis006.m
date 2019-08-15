function Large_sphere_analysis005(obs,ciedata2,xbar2,ybar2,zbar2,lambdaCie2,HighL,LowL,repstart)
%005 corrects the Z tristimulus values rather than the xy chromaticities of
%the averages, multiplicitively, hopefully allowing a more full analysis to
%take place

%%006 Created so that some pruning could take place losslessly in the
%%colour attributes of the uv plots. Will be making this section more
%%able to handle all observers data.


clc, clear, %close all

obs='DG';   %'LM', 'TR', or 'DG'

%% Load CIE data
ciefile = fullfile('C:','Users','ucesars','Dropbox','UCL','Data',...
    'Colour Standards','CIE colorimetric data','CIE_colorimetric_tables.xls');

if exist('xbar2','var')~=1
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
% figure('units','normalized','outerposition',[0 0 1 1])
% hold on
% axis('equal')
% title(obs)
% view(2)
% LowL=30;    LL=18-LowL/5;
% HighL=65;   HL=18-HighL/5;
%
% xlim([-30 30])
% ylim([-30 30])
% zlim([LowL-5 HighL+5])
%
% for i=1:length(files)
%     %ignore un-corrected DG 480 data
%     if strcmp(files(i).name,'480nm - time v2_3.mat')
%     else
%
%         LABdata=squeeze(mean(files(i).dataLABcal(:,HL:LL,:),3));
%         s=plot3(LABdata(2,:),LABdata(3,:),LABdata(1,:));
%
%         WL=(str2num(files(i).name(1:3)));
%         LCol=xyz2rgb(ciedata2(ciedata2(:,1)==WL,2:4));
%         LCol(LCol<0)=0;LCol(LCol>1)=1;
%         s.Color=LCol;
%
%         %drawnow
%         %pause(.25)
%         xlabel('a*')
%         ylabel('b*')
%         zlabel('L*')
%
%         %LABdata=reshape(files(i).dataLABcal(:,HL:LL,:),3,(LL-HL+1)*N);
%         %s=scatter3(LABdata(2,:),LABdata(3,:),LABdata(1,:),'r.');
%         %s.CData=[(str2num(files(i).name(1:3))-400)/300,0,0];
%         %s.SizeData=100;
%         %text(LABdata(2),LABdata(3),files(i).name)
%     end
% end
%
% % plot zero lines
% currentaxes=gca;
% plot([currentaxes.XLim],[0,0],'Color',[.8,.8,.8]);
% plot([0,0],[currentaxes.YLim],'Color',[.8,.8,.8]);
%% Plot all data in u'v'

LowL=35;    
HighL=60;   
repstart=3;

spectralambda=400:20:700;
LL=18-LowL/5;
HL=18-HighL/5;

fig=figure('DefaultAxesFontSize',24,'units','normalized','outerposition',[0 0 1 1],...
    'Color',[1 1 1]);
ax=axes();
hold on
axis('equal')
%title(sprintf('Obs-%s,L-%d:%d,repstart-%d',obs,LowL,HighL,repstart))
title({'Session-wise means of achromatic points, Observer #3',...
         'L* 35:60, Repititions 3:10'})

ciedata2_10001=interp1(lambdaCie2,ciedata2(:,2:4),380:.04:780,'spline');
xbar2_10001=interp1(lambdaCie2,xbar2,380:.04:780,'spline');
ybar2_10001=interp1(lambdaCie2,ybar2,380:.04:780,'spline');
zbar2_10001=interp1(lambdaCie2,zbar2,380:.04:780,'spline');
ubar=4*xbar2_10001./(xbar2_10001 + 15*ybar2_10001 + 3*zbar2_10001);
vbar=9*ybar2_10001./(xbar2_10001 + 15*ybar2_10001 + 3*zbar2_10001);
locus=scatter(ax,ubar,vbar,100,'filled');
LCol=xyz2rgb(ciedata2_10001);
LCol(LCol<0)=0;LCol(LCol>1)=1;
locus.CData=LCol;

% Scatter achromatic setting points (#1/2)
for j=4:13%1:length(spectralambda)
    for i=1:length(files)
        %takes the first data set for each wavelength, ignores repeats
        if strcmp(files(i).name(1:3),num2str(spectralambda(j)))
            s(j)=scatter(...
                mean(reshape(files(i).datauv(1,HL:LL,repstart:end),...
                (LL-HL+1)*(N-repstart+1),1)),...
                mean(reshape(files(i).datauv(2,HL:LL,repstart:end),...
                (LL-HL+1)*(N-repstart+1),1)),...
                'k','filled');
            break
        end
    end
end
xlabel('u''')
ylabel('v''')

% PR650 Sphere Measurements
load('Illumination in sphere');

% Calculate chromaticities and plot internal sphere measurements(400:20:700)
% (These are 2deg, which is appropriate for comparison but arguably
% incorrect generally for looking at these wide adapting fields)
for i=2:17
    
    X(i)=xbar2_101*spectra(:,i)*1000;
    Y(i)=ybar2_101*spectra(:,i)*1000;
    Z(i)=zbar2_101*spectra(:,i)*1000;
    
    u_prime(i)= 4*X(i)/(X(i) + 15*Y(i) + 3*Z(i));
    v_prime(i)= 9*Y(i)/(X(i) + 15*Y(i) + 3*Z(i));
    
end

scatter(u_prime(5:14),v_prime(5:14),'k','filled');

for i=5:14
    text(u_prime(i)+0.01,v_prime(i),num2str(spectralambda(i-1)),'FontSize', 15)
end


% Join achromatic points with corresponding sphere measurements
for j=4:13 %1:16
    plot([u_prime(j+1),s(j).XData] ,[v_prime(j+1),s(j).YData],'k');
end

axis('equal')
xlim([0 .7])
ylim([0.13 .6])
ax.YTick=[0.2:.1:.6];

% Make a subplot with a zoom of the achromatic points
axes('Position',[.56 .15 .25 .4])
box on
hold on
for j=4:13%1:length(spectralambda)
    for i=1:length(files)
        %takes the first data set for each wavelength, ignores repeats
        if strcmp(files(i).name(1:3),num2str(spectralambda(j)))
            s2(j)=scatter(...
                mean(reshape(files(i).datauv(1,HL:LL,repstart:end),...
                (LL-HL+1)*(N-repstart+1),1)),...
                mean(reshape(files(i).datauv(2,HL:LL,repstart:end),...
                (LL-HL+1)*(N-repstart+1),1)),...
                'k','filled');
            % text(0.2,0.5,'test');
            text(...
                mean(reshape(files(i).datauv(1,HL:LL,repstart:end),...
                (LL-HL+1)*(N-repstart+1),1))+0.0005,...
                mean(reshape(files(i).datauv(2,HL:LL,repstart:end),...
                (LL-HL+1)*(N-repstart+1),1))+0.0005,...
                files(i).name(1:3));
            if j>4 && j<14
                plot([s2(j-1).XData,s2(j).XData],[s2(j-1).YData,s2(j).YData],'k:')
            end
            break
        end
    end
end
axis('equal')
set(gca,'fontsize',10)
subx=xlim();
suby=ylim();

plot(ax,...
    [subx(1) subx(1) subx(2) subx(2) subx(1)],...
    [suby(1) suby(2) suby(2) suby(1) suby(1)],'k:');

plot_filename=fullfile('C:','Users','ucesars','Dropbox','UCL','Ongoing Work',...
    'Large Sphere','Figure 3');
img = getframe(gcf);
imwrite(img.cdata, [plot_filename, '.tiff']);

 %% Plot individual sessions in u'v'
% 
% % spectralambda=400:20:700;
% % LowL=35;
% % LL=18-LowL/5;
% % HighL=60;
% % HL=18-HighL/5;
% % repstart=3;
% 
% spectralambda=400:20:700;
% LowL=10;
% LL=18-LowL/5;
% HighL=85;
% HL=18-HighL/5;
% repstart=1;
% 
% fig=figure('DefaultAxesFontSize',24,'units','normalized','outerposition',[0 0 1 1],...
%     'Color',[1 1 1]);
% ax=axes();
% hold on
% axis('equal')
% %title(sprintf('Obs-%s,L-%d:%d,repstart-%d',obs,LowL,HighL,repstart))
% 
% %Adds spectral locus
% %-%
% ciedata2_10001=interp1(lambdaCie2,ciedata2(:,2:4),380:.04:780,'spline');
% xbar2_10001=interp1(lambdaCie2,xbar2,380:.04:780,'spline');
% ybar2_10001=interp1(lambdaCie2,ybar2,380:.04:780,'spline');
% zbar2_10001=interp1(lambdaCie2,zbar2,380:.04:780,'spline');
% ubar=4*xbar2_10001./(xbar2_10001 + 15*ybar2_10001 + 3*zbar2_10001);
% vbar=9*ybar2_10001./(xbar2_10001 + 15*ybar2_10001 + 3*zbar2_10001);
% locus=scatter(ax,ubar,vbar,100,'filled');
% LCol=xyz2rgb(ciedata2_10001);
% LCol(LCol<0)=0;LCol(LCol>1)=1;
% locus.CData=LCol;
% %-%
% 
% % Scatter achromatic setting points (#1/2)
% for i=5%1:length(files)%[9,18]%
%     s(i)=scatter(...
%         (reshape(files(i).datauv(1,HL:LL,repstart:end),...
%         (LL-HL+1)*(N-repstart+1),1)),...
%         (reshape(files(i).datauv(2,HL:LL,repstart:end),...
%         (LL-HL+1)*(N-repstart+1),1)),[],...
%          [linspace(.1,.9,(LL-HL+1)*(N-repstart+1))',...
%          linspace(.9,.1,(LL-HL+1)*(N-repstart+1))',...
%          zeros((LL-HL+1)*(N-repstart+1),1)],...
%          'filled');
%     title({'480nm achromatic points, Observer #3',...
%         'L* 10:85, Repititions 1:10'})
%     drawnow
%     %pause(3)
% 
% end
% 
% xlabel('u''')
% ylabel('v''')
% axis('equal')
% xlim([0 .7])
% ylim([0.13 .6])
% ax.YTick=[0.2:.1:.6];
% 
% % Make a subplot with a zoom of the achromatic points
% axes('Position',[.56 .15 .25 .4])
% box on
% hold on
% 
% 
% for i=5%1:length(files)%[9,18]%
%     s2(i)=scatter(...
%         (reshape(files(i).datauv(1,HL:LL,repstart:end),...
%         (LL-HL+1)*(N-repstart+1),1)),...
%         (reshape(files(i).datauv(2,HL:LL,repstart:end),...
%         (LL-HL+1)*(N-repstart+1),1)),[],...
%          [linspace(.1,.9,(LL-HL+1)*(N-repstart+1))',...
%          linspace(.9,.1,(LL-HL+1)*(N-repstart+1))',...
%          zeros((LL-HL+1)*(N-repstart+1),1)],...
%         'filled');
% end
% 
% axis('equal')
% set(gca,'fontsize',10)
% subx=xlim();
% suby=ylim();
% 
% plot(ax,...
%     [subx(1) subx(1) subx(2) subx(2) subx(1)],...
%     [suby(1) suby(2) suby(2) suby(1) suby(1)],'k:');
% 
% plot_filename=fullfile('C:','Users','ucesars','Dropbox','UCL','Ongoing Work',...
%     'Large Sphere','480_1');
% img = getframe(gcf);
% %imwrite(img.cdata, [plot_filename, '.tiff']);
% 
% % Plot individual sessions in u'v' (with data limitations
% 
% spectralambda=400:20:700;
% LowL=35;
% LL=18-LowL/5;
% HighL=60;
% HL=18-HighL/5;
% repstart=3;
% 
% fig=figure('DefaultAxesFontSize',24,'units','normalized','outerposition',[0 0 1 1],...
%     'Color',[1 1 1]);
% ax=axes();
% hold on
% axis('equal')
% %title(sprintf('Obs-%s,L-%d:%d,repstart-%d',obs,LowL,HighL,repstart))
% 
% %Adds spectral locus
% %-%
% ciedata2_10001=interp1(lambdaCie2,ciedata2(:,2:4),380:.04:780,'spline');
% xbar2_10001=interp1(lambdaCie2,xbar2,380:.04:780,'spline');
% ybar2_10001=interp1(lambdaCie2,ybar2,380:.04:780,'spline');
% zbar2_10001=interp1(lambdaCie2,zbar2,380:.04:780,'spline');
% ubar=4*xbar2_10001./(xbar2_10001 + 15*ybar2_10001 + 3*zbar2_10001);
% vbar=9*ybar2_10001./(xbar2_10001 + 15*ybar2_10001 + 3*zbar2_10001);
% locus=scatter(ax,ubar,vbar,100,'filled');
% LCol=xyz2rgb(ciedata2_10001);
% LCol(LCol<0)=0;LCol(LCol>1)=1;
% locus.CData=LCol;
% %-%
% 
% % Scatter achromatic setting points (#1/2)
% for i=5%1:length(files)%[9,18]%
%     s(i)=scatter(...
%         (reshape(files(i).datauv(1,HL:LL,repstart:end),...
%         (LL-HL+1)*(N-repstart+1),1)),...
%         (reshape(files(i).datauv(2,HL:LL,repstart:end),...
%         (LL-HL+1)*(N-repstart+1),1)),[],...
%          [linspace(.1,.9,(LL-HL+1)*(N-repstart+1))',...
%          linspace(.9,.1,(LL-HL+1)*(N-repstart+1))',...
%          zeros((LL-HL+1)*(N-repstart+1),1)],...
%          'filled');
%     title({'480nm achromatic points, Observer #3',...
%         'L* 35:60, Repititions 3:10'})
%     drawnow
%     %pause(3)
% 
% end
% 
% xlabel('u''')
% ylabel('v''')
% axis('equal')
% xlim([0 .7])
% ylim([0.13 .6])
% ax.YTick=[0.2:.1:.6];
% 
% % Make a subplot with a zoom of the achromatic points
% axes('Position',[.56 .15 .25 .4])
% box on
% hold on
% 
% 
% for i=5%1:length(files)%[9,18]%
%     s2(i)=scatter(...
%         (reshape(files(i).datauv(1,HL:LL,repstart:end),...
%         (LL-HL+1)*(N-repstart+1),1)),...
%         (reshape(files(i).datauv(2,HL:LL,repstart:end),...
%         (LL-HL+1)*(N-repstart+1),1)),[],...
%          [linspace(.1,.9,(LL-HL+1)*(N-repstart+1))',...
%          linspace(.9,.1,(LL-HL+1)*(N-repstart+1))',...
%          zeros((LL-HL+1)*(N-repstart+1),1)],...
%          'filled');
% end
% 
% axis('equal')
% set(gca,'fontsize',10)
% subx=xlim();
% suby=ylim();
% 
% plot(ax,...
%     [subx(1) subx(1) subx(2) subx(2) subx(1)],...
%     [suby(1) suby(2) suby(2) suby(1) suby(1)],'k:');
% 
% plot_filename=fullfile('C:','Users','ucesars','Dropbox','UCL','Ongoing Work',...
%     'Large Sphere','480_2');
% img = getframe(gcf);
% %imwrite(img.cdata, [plot_filename, '.tiff']);

 %% Hypothetical data
% 
% fig=figure('DefaultAxesFontSize',24,'units','normalized','outerposition',[0 0 1 1],...
%     'Color',[1 1 1]);
% ax=axes();
% hold on
% axis('equal')
% title({'Hypothetical ''Perfect'' Data for a ''trichromatic input''',...
%          'Surround Chromaticities'})
% 
% spectralambda=400:20:700;
% 
% %Adds spectral locus
% %-%
% ciedata2_10001=interp1(lambdaCie2,ciedata2(:,2:4),380:.04:780,'spline');
% xbar2_10001=interp1(lambdaCie2,xbar2,380:.04:780,'spline');
% ybar2_10001=interp1(lambdaCie2,ybar2,380:.04:780,'spline');
% zbar2_10001=interp1(lambdaCie2,zbar2,380:.04:780,'spline');
% ubar=4*xbar2_10001./(xbar2_10001 + 15*ybar2_10001 + 3*zbar2_10001);
% vbar=9*ybar2_10001./(xbar2_10001 + 15*ybar2_10001 + 3*zbar2_10001);
% locus=scatter(ax,ubar,vbar,100,'filled');
% LCol=xyz2rgb(ciedata2_10001);
% LCol(LCol<0)=0;LCol(LCol>1)=1;
% locus.CData=LCol;
% %-%
% 
% axis('equal')
% xlim([0 .7])
% ylim([0.13 .6])
% ax.YTick=[0.2:.1:.6];
% xlabel('u''')
% ylabel('v''')
% 
% drawnow;pause(1);
% plot_filename=fullfile('C:','Users','ucesars','Dropbox','UCL','Ongoing Work',...
%     'Large Sphere','Locus');
% img = getframe(gcf);
% imwrite(img.cdata, [plot_filename, '.tiff']);
% 
% % PR650 Sphere Measurements
% load('Illumination in sphere');
% 
% % Calculate chromaticities and plot internal sphere measurements(400:20:700)
% % (These are 2deg, which is appropriate for comparison but arguably
% % incorrect generally for looking at these wide adapting fields)
% for i=5:14
%     
%     X(i)=xbar2_101*spectra(:,i);
%     Y(i)=ybar2_101*spectra(:,i);
%     Z(i)=zbar2_101*spectra(:,i);
%     
%     u_prime(i)= 4*X(i)/(X(i) + 15*Y(i) + 3*Z(i));
%     v_prime(i)= 9*Y(i)/(X(i) + 15*Y(i) + 3*Z(i));
%     
%     scatter(u_prime(i),v_prime(i),'k','filled');
%     text(u_prime(i)+0.01,v_prime(i),num2str(spectralambda(i-1)),'FontSize', 15)
% end
% 
% axis('equal')
% xlim([0 .7])
% ylim([0.13 .6])
% ax.YTick=[0.2:.1:.6];
% 
% 
% drawnow;pause(1);
% plot_filename=fullfile('C:','Users','ucesars','Dropbox','UCL','Ongoing Work',...
%     'Large Sphere','Surround measurements');
% img = getframe(gcf);
% imwrite(img.cdata, [plot_filename, '.tiff']);
% 
% axis('equal')
% xlim([0 .7])
% ylim([0.13 .6])
% ax.YTick=[0.2:.1:.6];
% 
% estimatedata=[u_prime;v_prime].*0.1+0.2;
% estimatedata(2,:)=estimatedata(2,:)+0.2;
% scatter(estimatedata(1,5:14),estimatedata(2,5:14),'k', 'filled');
% 
% drawnow;pause(1);
% plot_filename=fullfile('C:','Users','ucesars','Dropbox','UCL','Ongoing Work',...
%     'Large Sphere','Achromatic Points');
% img = getframe(gcf);
% imwrite(img.cdata, [plot_filename, '.tiff']);
% 
% for j=5:14
%     plot([u_prime(j),estimatedata(1,j)] ,[v_prime(j),estimatedata(2,j)],'k');
% end
% 
% drawnow;pause(1);
% plot_filename=fullfile('C:','Users','ucesars','Dropbox','UCL','Ongoing Work',...
%     'Large Sphere','Points Joined');
% img = getframe(gcf);
% %imwrite(img.cdata, [plot_filename, '.tiff']);
end

