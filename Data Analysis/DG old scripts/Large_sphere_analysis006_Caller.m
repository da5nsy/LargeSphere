%clc, clear,
close all

obs1={'LM','TR','DG'};

repstart=3;
HighL=60;
LowL=35;
for obs2=1:3
    %Large_sphere_analysis006('DG',xbar2,ybar2,zbar2,lambdaCie2,HighL,LowL);
    %Large_sphere_analysis006(obs1{obs2},ciedata2,xbar2,ybar2,zbar2,lambdaCie2,HighL,LowL,repstart);
    Large_sphere_analysis007(obs1{obs2},ciedata2,xbar2,ybar2,zbar2,lambdaCie2,HighL,LowL,repstart);

end

%%
Large_sphere_analysis006('DG',ciedata2,xbar2,ybar2,zbar2,lambdaCie2,HighL,LowL,repstart)
