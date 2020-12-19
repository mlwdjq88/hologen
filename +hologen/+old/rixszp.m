N = 500;
idx = linspace(-1, 1, N);


[x, y] = meshgrid(idx);
z = 100;
lambda = .0002;
dx = .1;
ph1 = exp(2i*pi/lambda * sqrt((x + dx).^2 + y.^2 + z.^2));

ph2 = exp(2i*pi/lambda * sqrt((x - dx).^2 + y.^2 + z.^2));

am = (abs(2 + ph1 + ph2).^2) ;
imagesc(am);
%%
imagesc((am > max(am(:))/2).*logical(pinhole(N)))
axis square
