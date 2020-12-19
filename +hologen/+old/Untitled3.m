[x,y]=meshgrid(linspace(-5,5,100));
I=sin(x.^4+y.^2);
I=im2bw(I);
I=int32(I(:));
[L,num]=Cbwlabel(I,int32(100),int32(100));