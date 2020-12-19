function [L,num]=Cbwlabel(I,M,N)
I0=reshape(I,M,N);
[L,num]=bwlabel(I0,4);
L=int32(L(:));
num=int32(num);