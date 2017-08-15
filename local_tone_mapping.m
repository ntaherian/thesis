%% clustering image patches 
A = double(hdrread('nave.hdr'));
cform = makecform('srgb2lab');
%A = A(1:200, 1:200, :);
%A = ans;
%image_patch = [200 200 3];
image_patch = size(A);
img_size = size(A);
A_lab = applycform(A,cform);
ab = double(A_lab(:,:,2:3));
A_reshape = reshape(ab, img_size(1) * img_size(2),2);
smooth_lab = ksdensity(A_reshape);
K = size(findpeaks(smooth_lab),1);
%K = 95;
max_iters = 10;
w = 1 : image_patch(1) : img_size(1) - image_patch(1) + 1;
h = 1 : image_patch(2) : img_size(2) - image_patch(2) + 1;
NumOfPatches = size(w,2) * size(h,2);
centroids = zeros(NumOfPatches,K,3);
idx = zeros(image_patch(1)*image_patch(2),NumOfPatches);
n = 0;
for i = 1 : image_patch(1) : img_size(1) - image_patch(1) + 1
    for j = 1 : image_patch(2) : img_size(2) - image_patch(2) + 1
        X_patch = A(i:i+image_patch(1)-1,j:j+image_patch(2)-1,:);
        n = n + 1;
        X = reshape(X_patch, image_patch(1) * image_patch(2),3);
        max_X = max(X(:));
        additive = max_X/2000000000000000000000000000;
        X = sum(X.^2, 2).^.5;
        %X = mean(X,2);
        %X = max(X, [], 2);
        X = log(X + additive);
                
        initial_centroids = kMeansInitCentroids(X, K);
        [centroids_patch, idx_patch] = runkMeans(X, initial_centroids, max_iters, false);
        %centroids_patch_brightness = max(centroids_patch, [], 2);
        %centroids_patch_brightness = sum(centroids_patch.^2, 2).^.5;
        
        %centroids_patch = exp(centroids_patch(:,1));
        centroids_patch_brightness = exp(centroids_patch - additive);
        %centroids_patch_brightness = centroids_patch;
        C = sort(centroids_patch_brightness);
        I = zeros(size(C));
        for k = 1:K
            ii = find(centroids_patch_brightness(k) == C);
            I(k) = ii;
        end
        centroids_adjustment = repmat((I - 0.5) / (K) ./ centroids_patch_brightness, 1, 3);
        centroids(n,:,:) = centroids_adjustment;
        idx(:,n) = idx_patch;
    end
end
%% Image Compression
n = 0;
image_recovered = zeros([NumOfPatches image_patch]);
for i = 1 : image_patch(1) : img_size(1) - image_patch(1) + 1
    for j = 1 : image_patch(2) : img_size(2) - image_patch(2) + 1
        X_patch = A(i:i+image_patch(1)-1,j:j+image_patch(2)-1,:);
        n = n + 1;
        centroids_patch = reshape(centroids(n,:,:),K,3); 
        X = reshape(X_patch, image_patch(1) * image_patch(2),3);        
        adjustment_coeff = centroids_patch(idx_patch,:);
        adjustment_coeff = reshape(adjustment_coeff, image_patch);
        max_adjustment_coeff = max((adjustment_coeff(:)));
        adjustment_coeff = adjustment_coeff/max_adjustment_coeff;
        %adjustment_coeff = bfilter2(adjustment_coeff,5,1/700*max(image_patch(1),image_patch(2)));
        adjustment_coeff = imgaussfilt(adjustment_coeff,1/100*max(image_patch(1),image_patch(2)));
        %adjustment_coeff = medfilt3(adjustment_coeff)* max_adjustment_coeff;
        adjustment_coeff = adjustment_coeff * max_adjustment_coeff;
        X_recovered = X_patch .* adjustment_coeff;
        X_recovered = reshape(X_recovered, image_patch(1), image_patch(2), 3);
        image_recovered(n,:,:,:) = X_recovered;
    end
end
%% Putting image together 
image_clustered = zeros(size(A));
n = 1;
for i = 1 : image_patch(1) : img_size(1) - image_patch(1) + 1
    for j = 1 : image_patch(2) : img_size(2) - image_patch(2) + 1
        image_clustered(i:i+image_patch(1)-1,j:j+image_patch(2)-1,:) = image_recovered(n, :, :, :);
        n = n + 1;
    end
end
subplot(2, 2, 1); imagesc(A); title('Original');
subplot(2, 2, 2); imagesc(image_clustered); title('Tone-Mapped');

A_brightness = sum(A.^2,3).^.5;
subplot(2, 2, 3); hist(A_brightness(:),256); title('Original');

B_brightness = sum(image_clustered.^2,3).^.5;
subplot(2, 2, 4); hist(B_brightness(:),256); title('Tone-Mapped');