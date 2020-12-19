%unit:mm
lambda=13.5e-6;
I0=imread('C:\Users\T470P\OneDrive\ндуб\Interferometric measurement of etch-depth in etched phase shift masks\Figure\CXRO.bmp');
I1=rgb2gray(I0);
I1=im2bw(I1,0.7);
% figure,imshow(I1,[])
% I1(I1~=255)=0;
I2=(double(I1)-1)*pi/2;
%carrier of object
thita=6;
f=3;
R=tan(thita/180*pi).*f/2;
%number of rings
m=(sqrt(R.^2+f.^2)-f)/lambda;
%CCD size
cs=5;
%Magnification
Ma=250;
%imaging distance
lc=f.*Ma;
%FOV
FOV=cs*f./lc;
%NA
NA=sin(atan(R./f));
for T=0.013:-0.02:0.013
%grating pitch
% T=0.4;
%shear distance
s=lc.*(tan(asin(sin(0/180*pi)+lambda./T))-tan(asin(sin(0/180*pi)-lambda./T)));

%cal fringe vs incident angle
% for i=1:length(thita)
% res(i)=calFringeVSAngle(lambda,f,lc,T,thita(i));
% end
% res
Fn=2/T/Ma*cs;

[x,y]=meshgrid(linspace(-2.5,2.5,500));
object=zeros(500);
% object(1:100,1:256)=pi;
% object(150:200,231:281)=0.5*pi;
% object(300:320,246:266)=0.25*pi;
% object(400:500,100:200)=0.4*pi;
% object(400:500,300:400)=0.4*pi;
%  object(236:276,236:276)=-pi;
%  obs= circshift(object,50,2)-circshift(object,-50,2);
 object=I2;
% y=y+offset;
num=1;
% figure(2),imagesc(object)
for i=0:0.05:0.5
th=i/180*pi;
ph0=asin(sin(0)+lambda/T);
ph1=asin(sin(th)+lambda/T);
ph2=asin(sin(th)-lambda/T);
offset0=lc*tan(ph0);
offsetpixel=round(offset0/cs*500);
fl=29*f;
offsetonHolo=f*tan(th);
offset1=fl*tan(ph1)+offsetonHolo;
offset2=fl*tan(ph2)+offsetonHolo;
object1=object;
object1(:,1:500-offsetpixel)=object1(:,1+offsetpixel:500);
object1(:,500-offsetpixel+1,end)=0;
object2=object;
object2(:,1+offsetpixel:500)=object2(:,1:500-offsetpixel);
object2(:,1:offsetpixel)=0;
dobj=object1-object2;
interferogram=1+cos(2*pi/lambda*(sqrt((x-fl*tan(ph1)).^2+(y+(lc-fl)*tan(thita/180*pi)).^2+...
    (lc-fl)^2)-sqrt((x-fl*tan(ph2)).^2+(y+(lc-fl)*tan(thita/180*pi)).^2+(lc-fl)^2)+lambda/T*offsetonHolo)+dobj);
fi=fftshift(fft2(interferogram));
fi2=zeros(500);
Ns=125;
fi2(251-Ns:251+Ns,251-Ns:251+Ns)=fi(251-Ns:251+Ns,251-Ns:251+Ns);
interferogram2=abs(ifft2(ifftshift(fi2)));
mask=zeros(500);
mask((x-offset0).^2+y.^2<(cs/2)^2)=-1;
mask((x+offset0).^2+y.^2<(cs/2)^2)=-1;
mask((x-offset0).^2+y.^2<(cs/2)^2&(x+offset0).^2+y.^2<(cs/2)^2)=1;
interferogram=interferogram2.*mask;
interferogram(interferogram<0)=1;
angs(num)=i;
s(num)=interferogram(188,170);
s2(num)=interferogram(170,170);
Is(num,:)=interferogram(mask==1);
num=num+1;
figure(1),imshow(interferogram,[])
drawnow;
end
end
figure(2),plot(angs,s),hold on
plot(angs,s2),hold on
mask2=mask;
mask2(mask~=1)=0;
mask2(mask==1)=255;
delta=linspace(0,2*pi,num-1);
[phase,deltas,k] = RandomShift(Is,delta,20);
ph=zeros(500);
ph(mask==1)=phase;
phs = UnwrapPhaseBySortingReliabilityWithMask(ph,mask2); 
phs(mask~=1)=NaN;
phs=DelTilt(phs);
figure(3),imagesc(x(1,:)*100/25,y(:,1)*100/25,phs);

%% fft
ft=fftshift(fft2(interferogram-mean(interferogram(:))));

win=hanning2d(ones(500), [250,250+106], 1.6*106);
ft2=ft.*win;
% ft2((x).^2+y.^2<0.1^2)=0;
% ft2((x-1.06).^2+y.^2>0.75^2)=0;
ph2=ifft2(fftshift(ft2));
ph3=atan2(imag(ph2),real(ph2));
ph4 = UnwrapPhaseBySortingReliabilityWithMask(ph3,mask2); 
ph4(mask2==0)=NaN;
ph4=DelTilt(ph4);
ph4=DelTilt(ph4);
figure(4),imagesc(x(1,:)*100/25,y(:,1)*100/25,-ph4);



