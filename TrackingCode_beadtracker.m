%%% Script that analyzes the displacement of beads
%%% Uses image files from Camera Control
%%% Arash Tajik, August 2014, Rev 1
%%% Fazlur Rashid, June 2024

clc;
clear;
close all;

pathname='C:\\';

[filename1,pathname]=uigetfile([pathname,'*.*'],'select first image in sequence');
[filename2,pathname]=uigetfile([pathname,'*.*'],'select last image in sequence');


cd(pathname);
mkdir(pathname,'AnalysisBEAD');

I=imread([pathname,filename1]);


[I2 rect]=imcrop(I);
imshow(I2);

k1 = strfind(filename1, '_');
k2 = strfind(filename1, '+');
k3= strfind(filename1, '-');

k=[k1 k2 k3];
k=sort(k);
k1=k(1);
k2=k(2);


CounterStart=str2num(filename1([k1+1:k2-1]));
CounterEnd=str2num(filename2([k1+1:k2-1]));


PositionList=[];
mkdir([pathname,'AnalysisBEAD'],'Crop');
mkdir([pathname,'AnalysisBEAD'],'Centroid');

cd([pathname,'AnalysisBEAD']);
writerObj = VideoWriter('Video.avi');
open(writerObj);


BeadCenter=[];
time=[];

for i=[CounterStart:CounterEnd]
    filename=['*',num2str(i),'*.*'];
    cd(pathname);
    filenameX=dir([filename]);          %% Find images in range
    
    if isempty(filenameX)
        continue
    end
    time=[time i];
    filename=filenameX.name;
    I=imread([pathname,filename]);
    I2=imcrop(I,rect);
    
    cd([pathname,'AnalysisBEAD','/Crop']);
    imwrite(I2,[filename1(1:k1),num2str(i),'_crop.tiff'],'tiff')

    [centers radii] = imfindcircles(I2,[1 35],'ObjectPolarity','dark');
    %[centers radii] = imfindcircles(I2,2,'ObjectPolarity','dark');


    xc=centers(1);
    yc=centers(2);
    
    BeadCenter=[BeadCenter;[xc yc]];
    
    
    figure(i);imshow(I2);hold on;plot(xc,yc,'ro');
    hold on; viscircles(centers, radii,'EdgeColor','b');
    cd([pathname,'AnalysisBEAD','/Centroid']);
    saveas(gcf,[filename1(1:k1),num2str(i),'_Centroid.fig']);
    saveas(gcf,[filename1(1:k1),num2str(i),'_Centroid.tiff']);
    frame = getframe;
    writeVideo(writerObj,frame);
    close(gcf)

end

close(writerObj);

figure; plot(BeadCenter(:,1),BeadCenter(:,2),'bo-')
cd([pathname,'AnalysisBEAD']);
saveas(gcf,'BeadDisplacement.tiff');
saveas(gcf,'BeadDisplacement.fig');

drift=BeadCenter(1,:)-BeadCenter(end,:);
drift=drift/(size(BeadCenter,1)-1);
BeadCenterD2=BeadCenter+[drift(1)*[0:(size(BeadCenter,1)-1)]' drift(2)*[0:(size(BeadCenter,1)-1)]']
figure; plot(BeadCenterD2(:,1),BeadCenterD2(:,2),'ro-')
cd([pathname,'AnalysisBEAD']);
saveas(gcf,'BeadDisplacementwoDrift.tiff');
saveas(gcf,'BeadDisplacementwoDrift.fig');


dxdy=BeadCenterD2(2:end,:)-BeadCenterD2(1:end-1,:);
disp=sqrt(dxdy(:,1).^2+dxdy(:,2).^2);
% ang=atan2(dxdy(:,2),dxdy(:,1));
ang=atan(dxdy(:,2)./dxdy(:,1));
figure;plot(disp)
cd([pathname,'AnalysisBEAD']);
saveas(gcf,'BeadDisplacementAmplitudewoDrift.tiff');
saveas(gcf,'BeadDisplacementAmplitudewoDrift.fig');
figure;rose(ang(disp>(mean(disp))),[0:15:360])
cd([pathname,'AnalysisBEAD']);
saveas(gcf,'BeadDisplacementAnglewoDrift.tiff');
saveas(gcf,'BeadDisplacementAnglewoDrift.fig');

pathname


cd([pathname,'AnalysisBEAD']);
save AllVariableBEAD.mat






   





        
        
        
    








