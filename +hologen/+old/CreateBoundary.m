function CreateBoundary(outputFile,xy,N)
[polyPre,polyForm,polyPost]=polyDef(N);
coords=[polyPre;polyForm;xy(:);polyPost];
fwrite(outputFile,coords, 'int32','b');


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