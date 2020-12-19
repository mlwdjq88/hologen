function error=CalError(x,y,delta,f,lambda,th,db)
lambda=lambda*db;
delta=delta*db;
f=f*db;
E=exp(1i*2*pi/lambda*sqrt((x-delta/2).^2+y.^2+f.^2))+...
            exp(1i*2*pi/lambda*sqrt((x+delta/2).^2+y.^2+f.^2));
I=abs(2+E);
I=I.^2;
error=abs(I-th);