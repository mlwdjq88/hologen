function d=Plustest(a)
N=length(a);
b=a;
c=a;
for i=1:N^2
        c(i)=a(i)+b(i);
end
d=c(1);