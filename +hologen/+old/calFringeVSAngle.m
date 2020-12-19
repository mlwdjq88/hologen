function res=calFringeVSAngle(lambda,f,z,T,th0)
%%parameter setting

th=th0/180*pi+linspace(-1/180*pi,1/180*pi,1000);%set incident angle range
% [th T]=meshgrid(th,T);
th1=asin(sin(-th)+lambda./T);
th2=asin(sin(-th)-lambda./T);

% Part 1
ph1=pi./cos(th);
% Part 2
delta=z*f.*(tan(th1)-tan(th2))/(z-f);
ph2=-2*pi.*delta./T;
% Part 3
delta0=0 ;%central offset
deltal1=f*(tan(th0/180*pi)-tan(th))+delta0;
deltal2=(f+lambda/2)*(tan(th0/180*pi)-tan(th))+delta0;
p1=atan(f./(delta/2+deltal1-f*tan(th1)));
p2=atan(f./(-f*tan(th2)-delta/2+deltal2));
l1=z./sin(p1);
l2=z./sin(p2);
ph3=2*pi*(l2-l1)/lambda;
% Part 4
ph4=2*pi/lambda*(sqrt((delta/2+deltal1).^2+f^2)-sqrt((delta/2-deltal2).^2+f^2));
% total phase differnce
% ph1=0;
ph=ph1+ph2+ph3+ph4;
% fringe intensity
Intensity=1+cos(ph);

figure(1),plot(th/pi*180,Intensity)

d0=z./tan(p1);
deltas=delta-lambda/2*tan(th);
d1=z./(f./(delta/2-f*tan(th2)+deltal1));
d2=z./(f./(-delta/2-f*tan(th1)+deltal2))+delta;
figure(2),plot(th/pi*180,[d1-d0])

res=(max(th)-min(th))/pi*180/((max(ph)-min(ph))/2/pi);%degree/period
spatialfringenum=30e-3/T;