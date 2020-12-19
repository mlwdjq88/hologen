% clc
% clear;
clear I3
I=double(uimread());
[M,N,~]=size(I);
num=6;
Is=zeros(M,N,6);
for i=1:6
    for j=1:6
    Is(:,:,i)=Is(:,:,i)+I(:,:,(i-1)*6+j);
    end
end
    
figure(1),imshow(Is(:,:,1),[]);
[x,y]=ginput(2);
close(1);
mask=zeros(M,N);
mask(y(1):y(2),x(1):x(2))=1;
% figure(2),imshow(mask,[])
I2=reshape(Is,M*N,num);
for i=1:num
    temp = I2(:,i);
    I3(i,:)= temp(mask(:)==1);
end
delta=linspace(-pi,pi-2*pi/num,num);
limit=30;
[phase,deltas,k] = RandomShift(I3,delta,limit);
deltas=unwrap(deltas)/pi*180;
p=polyfit(1:num,deltas,1);
ph1=zeros(M,N);
if p(1)<0
    ph1(mask==1)=-phase;
else
    ph1(mask==1)=phase;
end

mask=double(mask)*255;
ph2 = UnwrapPhaseBySortingReliabilityWithMask(ph1,mask);  %相位解包
% ph2=ph1;
ph2(mask==0)=NaN;
% Iout=ph2/2/pi;
% mask2=zeros(M,N);
% mask2(y(3):y(4),x(3):x(4))=1;
% % ph3=ph2;
% % ph3(mask2~=1)=NaN;
% ph3=DelTilts(ph2,mask2);
phs=ph2;
try
 dp=phs-phs2;dp=dp-mean(dp(~isnan(dp)));
 dp(isnan(dp))=0;
fdp=fftshift(fft2(dp));
[M,N]=size(fdp);
fdps=zeros(M,N);
r=10;
fdps(M/2-r:M/2+r,N/2-r:N/2+r)=fdp(M/2-r:M/2+r,N/2-r:N/2+r);
%  figure(5),mesh(dp/2/pi)
ph4=ifft2(ifftshift(fdps));
figure(5),mesh(real(ph4/2/pi))
catch
end
% figure(2),mesh(ph2);