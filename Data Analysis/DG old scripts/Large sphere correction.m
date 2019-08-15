%plot XYZ of DG large sphere data for the 4 480s
%See if there's a simple translation

%Divide Z/1.7 and we seem to be on to a winner.

%Do first sections of large_sphere_analysis004, up to 'CORRECTION'

figure, hold on
for i=5:8%1:length(files)
    XYZdata=reshape(files(i).dataXYZ,3,160);
    if i==7
        s=scatter3(XYZdata(1,:),XYZdata(3,:)/1.7,XYZdata(2,:),'r');
    else
        s=scatter3(XYZdata(1,:),XYZdata(3,:),XYZdata(2,:),'filled');
    end
end

xlabel('X');
ylabel('Z');
zlabel('Y');