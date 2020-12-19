function I=getIntensity(xs,ys,delta,f,lambda)
E=exp(1i*2*pi/lambda*sqrt((xs-delta/2).^2+ys.^2+f.^2))+...
    exp(1i*2*pi/lambda*sqrt((xs+delta/2).^2+ys.^2+f.^2));
A=ones(size(E));
I=A.*abs(2+E);
I=I.^2;