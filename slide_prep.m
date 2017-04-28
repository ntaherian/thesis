%% Slide 1
% Try log of intensity.
% Use X and Y as a feature for K-Means.
imshow(A/max(A(:)), [])

Z = max(X_patch, [], 3);
%Z = reshape(X, image_patch(1:2))
figure;
imshow(Z/max(Z(:)), [])
figure;
hist(X/max(X), 50);

%% Slide 2 - 1
Z = reshape(idx_patch, image_patch(1:2));
imagesc(Z)

%% Slide 2 - 1
X2 = max(X, [], 2);
X2 = X2 / max(X2);
X4 = [];
x_axis = linspace(0,1,200);
for k = 1:K
    X3 = X2(idx_patch == k);
    X3 = hist(X3, x_axis);
    X4 = [X4; X3];
end
bar(x_axis, X4', 7)
%% Slide 3
X2 = max(X_recovered, [], 3);
X2 = X2(:);
%X2 = X2 / max(X2);
X4 = [];
x_axis = linspace(0,2,200);
for k = 1:K
    X3 = X2(idx_patch == k);
    X3 = hist(X3, x_axis);
    X4 = [X4; X3];
end
bar(x_axis, X4', 7)
hold all;
for k = 1:K
    plot([k/K, k/K], [-100, 1e5])
end

%% Slide 4
subplot(2,2,1)
imagesc(adjustment_coeff(:,:,1))
title('Original Coefficients')

subplot(2,2,2)
imagesc(imgaussfilt(adjustment_coeff(:,:,1), context.sigma))
%imagesc(bfilter2(adjustment_coeff(:,:,1), 7,2))
title('Smoother Coefficients')

subplot(2,2,3)
imagesc(bfilter2(adjustment_coeff(:,:,1), 7,2))
