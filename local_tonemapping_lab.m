%% clustering image patches 
tic
A = double(hdrread('2.hdr'));
%A = A(1:200, 1:200, :);
%image_patch = [200 200 3];
image_patch = size(A);
img_size = size(A);
cform = makecform('srgb2lab');
%K = 89;
max_iters = 6;
w = 1 : image_patch(1) : img_size(1) - image_patch(1) + 1;
h = 1 : image_patch(2) : img_size(2) - image_patch(2) + 1;
NumOfPatches = size(w,2) * size(h,2);
%centroids = zeros(NumOfPatches,K,3);
idx = zeros(image_patch(1)*image_patch(2),NumOfPatches);
n = 0;
for i = 1 : image_patch(1) : img_size(1) - image_patch(1) + 1
    for j = 1 : image_patch(2) : img_size(2) - image_patch(2) + 1
        X_patch = A(i:i+image_patch(1)-1,j:j+image_patch(2)-1,:);
        lab_X_patch = applycform(X_patch,cform);
        ab = double(lab_X_patch(:,:,2:3));
        n = n + 1;
        X_K = reshape(ab, image_patch(1) * image_patch(2),2);
        X = reshape(X_patch, image_patch(1) * image_patch(2),3);
        smooth_lab_hist = ksdensity(X_K);
        K = size( findpeaks(smooth_lab_hist),1);
%         hsv = rgb2hsv(X);
%         smooth_hue_hist = ksdensity(hsv(:,1));
%         K = size(findpeaks(smooth_hue_hist),2);
        centroids = zeros(NumOfPatches,K,3);
        %X = sum(X.^2, 2).^.5;
        X = mean(X,2);
        %X = max(X, [], 2);
        %X = log(X);
        
%         [cor1, cor2] = meshgrid(1:img_size(1),1:img_size(2));
%         cor1 = cor1/img_size(1);
%         cor1 = cor1(:);
%         cor2 = cor2/img_size(2);
%         cor2 = cor2(:);
%         X = [X, cor1, cor2];
                
        initial_centroids = kMeansInitCentroids(X, K);
        [centroids_patch, idx_patch] = runkMeans(X, initial_centroids, max_iters, false,1);
        %centroids_patch_brightness = max(centroids_patch, [], 2);
        %centroids_patch_brightness = sum(centroids_patch.^2, 2).^.5;
        if (any(isnan(centroids_patch)) == true)
            idx_nan = isnan(centroids_patch);
            b = centroids_patch(~isnan(centroids_patch));
            centroids_patch(idx_nan) = mean(b);
        end
        %centroids_patch = exp(centroids_patch(:,1));
        centroids_patch_brightness = centroids_patch;
        C = sort(centroids_patch_brightness);
        I = zeros(size(C));
        for k = 1:K
            ii = find(centroids_patch_brightness(k) == C);
            if (length(ii) ~= 1)
               I(k) = ii(1);
            else
               I(k) = ii;
            end
        end
        centroids_adjustment = repmat((I - 0.5) / (K) ./ centroids_patch_brightness, 1, 3);
        %centroids_adjustment = histeq(centroids_adjustment,K);
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
        idx_patch = idx(:, n);
        adjustment_coeff = centroids_patch(idx_patch,:);
        adjustment_coeff = reshape(adjustment_coeff, image_patch);
        %adjustment_coeff = medfilt3(adjustment_coeff);
        max_adjustment_coeff = max(adjustment_coeff(:));
        adjustment_coeff = adjustment_coeff/max_adjustment_coeff;
        %adjustment_coeff = wiener2(adjustment_coeff);
        adjustment_coeff = bfilter2(adjustment_coeff,3, 1/50*max(image_patch(1),image_patch(2)));
        %adjustment_coeff = imgaussfilt3(adjustment_coeff,1/50*max(image_patch(1),image_patch(2)));
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
toc