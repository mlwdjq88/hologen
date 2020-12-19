function [polyPre,polyForm,polyPost]=polyDef()
polyPre   = [0, 4, 8, 0, ...
    0, 6, 13, 2, 0, 1, ...
    0, 6, 14, 2, 0, 0];
polyForm  = [0, 44, 16, 3];
polyPost    = [0, 4, 17, 0];
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