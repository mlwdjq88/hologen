function GenGDSfeature(B)
filename='features.gds';
filenamelib='noname';
filenamestr='features';
outputFile=fopen(filename,'wb');
HEADER=[0, 6, 0, 2, 0, 7];
BGNLIB = [ 0, 28, 1, 2, 230, 43, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0,...
    230, 43, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0];
Llib=length(filenamelib)+4;
Llib=dec2hex(Llib,4);
L(1)=hex2dec(Llib(1:2));
L(2)=hex2dec(Llib(3:4));
LIBNAME=[L,2,6,double(filenamelib)];
UNITS=[0, 20, 3, 5, 61, 104, 219, 139, 172, 113, 12, 180, 56, 109,...
    243, 127, 103, 94, 246, 236];
BGNSTR=[0, 28, 5, 2, 0, 114, 0, 4, 0, 17, 0, 13,...
    0, 22, 0, 56, 0, 114, 0, 4, 0, 17, 0, 13, 0, 22, 0, 56];
Lstr=length(filenamestr)+4;
Lstr=dec2hex(Lstr,4);
L(1)=hex2dec(Lstr(1:2));
L(2)=hex2dec(Lstr(3:4));
STRNAME=[L,6,6,double(filenamestr)];
gdspreamble=[HEADER,BGNLIB,LIBNAME,UNITS,BGNSTR,STRNAME];
fwrite(outputFile,gdspreamble, 'uint8' );
for p=1:10
    for q=1:10
        for i=1:29
            N=size(B{i},1);
            % N=5;%多边形顶点数
            xy=[B{i};B{i}(1,:)];
            xyrot90=xy;
            xyrot90(:,1)= -xy(:,2);
            xyrot90(:,2)= xy(:,1);
            xys=xyrot90';
            if p<=5
                xys(1,:)=xys(1,:)+(p-1)*80000000;
            else
                xys(1,:)=xys(1,:)+(p-1)*80000000+20000000;
            end
            xys(2,:)=xys(2,:)+(q-1)*80000000;
            [polyPre,polyForm,polyPost]=polyDef(N);
            coords=[polyPre;polyForm;xys(:);polyPost];
            fwrite(outputFile,coords, 'int32','b');
        end
    end
end
gdsPost = [0, 4, 7, 0, 0, 4, 4, 0];
fwrite(outputFile,gdsPost , 'uint8' );
fclose(outputFile);
winopen(filename);


function [polyPre,polyForm,polyPost]=polyDef(N)
polyPre   = [0, 4, 8, 0, 0, 6, 13, 2, 0, 1, 0, 6, 14, 2, 0, 0];
polyPost    = [0, 4, 17, 0];
Llib=4+8*(N+1);
Llib=dec2hex(Llib,4);
L(1)=hex2dec(Llib(1:2));
L(2)=hex2dec(Llib(3:4));
polyForm  = [L, 16, 3];
str=dec2hex(polyPre,2);
str=str';
str=str(:);
polyPres(1,:)=str(1:8)';
polyPres(2,:)=str(9:16)';
polyPres(3,:)=str(17:24)';
polyPres(4,:)=str(25:32)';
polyPre=hex2dec(polyPres);
str=dec2hex(polyForm);
str=str';
str=str(:);
polyForm=str';
polyForm=hex2dec(polyForm);
str=dec2hex(polyPost);
str=str';
str=str(:);
polyPost=str';
polyPost=hex2dec(polyPost);