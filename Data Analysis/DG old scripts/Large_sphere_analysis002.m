%% Loads data

clc, clear, %close all

rootdir = fullfile('C:','Users','ucesars','Dropbox','UCL','Data','Large Sphere','Results - Oct 2016');
%rootdir = fullfile('C:','Users','ucesars','Dropbox','UCL','Data','Large Sphere','Lindsay time series - Sep 2012');
%rootdir = fullfile('C:','Users','ucesars','Dropbox','UCL','Data','Large Sphere','Tania  time series - Apr 2013');

cd(rootdir)
files= dir('*.mat');

N = 10;                             % number of repetitions over time
LN = 16;                            % number of lightness levels per repeat

for j=1:length(files)
    
    load(fullfile(rootdir,files(j).name));  % load experimental results
    %LABmatchFull(:,:,:,j)=LABmatch;         % pull in to composite
    files(j).data=LABmatch;     
    files(j).data=RGBmatch;
    
end

clear LABmatch j

%% Plot average of runs end-3 and end-1

figure,
axis([-20 25 -45 15])
axis('equal')

highL=70;   %max 85
lowL=20;    %min 10
scaler=2;   %size of points

hold on
for i=1:length(files)
    for j=18-highL/5:18-lowL/5       

        A(i,j)=mean(files(i).data(2,j,end-3:end-1));
        B(i,j)=mean(files(i).data(3,j,end-3:end-1));
        if strcmp(files(i).date(1:6),'22-Novb')|...
                strcmp(files(i).date(1:6),'25-Novb')
            sc1= scatter(A(i,j),B(i,j),(17-j)^scaler,'b','filled');
            pause(0.05)
        else           
            sc1= scatter(A(i,j),B(i,j),(17-j)^scaler,'r','filled');
            pause(0.05)
            sc1= scatter(A(i,j),B(i,j),(17-j)^scaler,([.6,.6,.7]),'filled');
        end
    end
    text(median(median(files(i).data(2,:,:)))-30, ...
        median(median(files(i).data(3,:,:))),...
        files(i).name)
    pause(0.4)
end

title(sprintf('Median across end-3 : end-1 sessions, for L=%d to L=%d',lowL,highL));
xlabel('A')
ylabel('B')

% plot zero lines
currentaxes=gca;
plot([currentaxes.XLim],[0,0],'Color',[.8,.8,.8]);
plot([0,0],[currentaxes.YLim],'Color',[.8,.8,.8]);