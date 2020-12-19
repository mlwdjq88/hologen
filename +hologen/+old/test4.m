clc
clear all
% mex CgetPolygonCoords.cpp
% mex CdownSampling.cpp
% tic,
%% parameter setting
db=10000000; %set unit to anstrom
lambda=13.5e-6;
delta=0.008;
f=3;
NA=0.0875;
R=f*tan(asin(NA));
Nx=500;
Ny=200;
dire=1;
AccuCtrl=0.001;
Mx=1;
My=1;
for p=1:Mx
    for q=1:My
        getCoords_old(lambda,delta,f,R,db,Nx,Ny,p,q,dire,AccuCtrl);
    end
end
% toc
%% generate GDS
tic
filename='test.gds';
outputFile=fopen(filename,'wb');
gdsPost=initGDS(outputFile);
for p=1:Mx
    for q=1:My
        filename=strcat('test',num2str(p),'&',num2str(q),'.cor');
        fid=fopen(filename,'rb');
        allcoords=fread(fid,'int32','l');
        fwrite(outputFile,allcoords, 'int32','b');
        fclose(fid);
        % error analysis
        Ex(:,1)=allcoords(6:16:end);
        Ex(:,2)=allcoords(8:16:end);
        Ex(:,3)=allcoords(10:16:end);
        Ex(:,4)=allcoords(12:16:end);
        Ex(:,5)=allcoords(14:16:end);
        Ey(:,1)=allcoords(7:16:end);
        Ey(:,2)=allcoords(9:16:end);
        Ey(:,3)=allcoords(11:16:end);
        Ey(:,4)=allcoords(13:16:end);
        Ey(:,5)=allcoords(15:16:end);
        Ex=Ex(:);
        Ey=Ey(:);
        error=CalError(Ex,Ey,delta,f,lambda,8,db);
        error(error>3*mean(error))=[];
    error(error>mean(error)+3*std(error))=[];
      er(p,q)=mean(error)/8;
         clear Ex Ey
%         figure(5),plot(error),hold on
    end
end
fwrite(outputFile,gdsPost,'uint8');
fclose(outputFile);
toc
delete('*.cor');
winopen('test.gds')