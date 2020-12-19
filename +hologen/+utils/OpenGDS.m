function outputFile=OpenGDS(filename,filenamestr)
filenamelib='noname';
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

