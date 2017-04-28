function rgb = yuv2rgb(yuv)
y = yuv(:,:,1);
u = yuv(:,:,2);
v = yuv(:,:,3);

r = 1 * y + 0 * u + 1.13983 * v;
g = 1 * y - 0.39465 * u - 0.58060 * v;
b = 1 * y + 2.03211 * u - 0 * v;

rgb = zeros(size(yuv));
rgb(:,:,1) = r;
rgb(:,:,2) = g;
rgb(:,:,3) = b;