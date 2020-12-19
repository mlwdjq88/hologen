function maskGDS(filename,x,y,offsetx,offsety)
outputFile=fopen(filename,'wb');
nameof6=double(filename(1:6));
gdspreamble = [ 0, 6, 0, 2, 0, 7, 0, 28, 1, 2, 230, 43, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0,...
    230, 43, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 10, 2, 6, 110, 111, 110, 97,...
    109, 101, 0, 20, 3, 5, 61, 104, 219, 139, 172, 113, 12, 180, 56, 109,...
    243, 127, 103, 94, 246, 236, 0, 28, 5, 2, 0, 114, 0, 4, 0, 17, 0, 13,...
    0, 22, 0, 56, 0, 114, 0, 4, 0, 17, 0, 13, 0, 22, 0, 56,...
    0, 10, 6, 6];
[polyPre,polyForm,polyPost]=polyDef();
fwrite(outputFile,gdspreamble, 'uint8' );
fwrite(outputFile,nameof6, 'uint8' );
[LC,~]=size(x);
coords=zeros(16,LC);
coords(1:5,:)=repmat([polyPre;polyForm],1,LC);
coords(6:2:14,:)=x'+offsetx;
coords(7:2:15,:)=y'+offsety;
coords(16,:)=repmat(polyPost,1,LC);
fwrite(outputFile,coords, 'int32','b');
gdsPost     = [0, 4, 7, 0, 0, 4, 4, 0];
fwrite(outputFile,gdsPost , 'uint8' );
fclose(outputFile);


function [polyPre,polyForm,polyPost]=polyDef()
polyPre   = [0, 4, 8, 0, 0, 6, 13, 2, 0, 1, 0, 6, 14, 2, 0, 0];
polyPost    = [0, 4, 17, 0];
polyForm  = [0, 44, 16, 3];
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