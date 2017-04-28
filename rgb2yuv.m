function yuv = rgb2yuv(rgb)
r = rgb(:,:,1);
g = rgb(:,:,2);
b = rgb(:,:,3);

y = 0.299 * r + 0.587 * g + 0.114 * b;
u = -0.14713 * r - 0.28886 * g + 0.436 * b;
v = 0.615 * r - 0.51499 * g - 0.10001 *b;

yuv = zeros(size(rgb));
yuv(:,:,1) = y;
yuv(:,:,2) = u;
yuv(:,:,3) = v;