function angle=calAngle(x0,y0,x1,y1,x2,y2)
A=[x1-x0,y1-y0];
B=[x2-x1,y2-y1];
y=det([A;B]);
% y=(x1-x0).*(y2-y1)-(x2-x1).*(y1-y0);
x=dot(A,B);
% x=(x1-x0).*(x2-x1)+(y1-y0).*(y2-y1);
angle=atan2(y,x);