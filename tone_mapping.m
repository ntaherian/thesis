function [out, context] = tone_mapping(A, context)

cform = makecform('srgb2lab');

if context.initialize 
    size_frame = size(A);
    A_lab = applycform(A,cform);
    ab = double(A_lab(:,:,2:3));
    A_reshape = reshape(ab, size_frame(1) * size_frame(2),2);
    smooth_lab = ksdensity(A_reshape);
    context.K = size(findpeaks(smooth_lab),1);
end
% if context.initialize 
%         size_frame = size(A);
%         A_reshape = reshape(A, size_frame(1) * size_frame(2),3);
%         max_A = max(A_reshape(:));
%         A_reshape = A_reshape/max_A;
%         hsv = rgb2hsv(A_reshape);
%         smooth_hue_hist = ksdensity(hsv(:,1));
%         context.K = size(findpeaks(smooth_hue_hist),2);
% end
K = context.K;
%K =45;
image_patch = size(A);
img_size = size(A);

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
        additive = max_X/50000;
%         hsv = rgb2hsv(X);
%         smooth_hue_hist = ksdensity(hsv(:,1));
%         K = size(findpeaks(smooth_hue_hist),2);
%         centroids = zeros(NumOfPatches,K,3);
        X = sum(X.^2, 2).^.5;
        %X = mean(X,2);
        %X = min(X, [], 2);
        %X = max(X, [], 2);
        X = log(X + additive);

        
        %         [cor1, cor2] = meshgrid(1:img_size(1),1:img_size(2));
        %         cor1 = cor1/img_size(1);
        %         cor1 = cor1(:);
        %         cor2 = cor2/img_size(2);
        %         cor2 = cor2(:);
        %         X = [X, cor1, cor2];
        
        if (context.initialize)
            context.initialize = false;
            max_iters = context.initial_iters;
            initial_centroids = kMeansInitCentroids(X, K);
            kmeans_alpha = 1;
        else
            max_iters = context.iters;
            initial_centroids = context.centroids;
            kmeans_alpha = 0.2;
        end

        [centroids_patch, idx_patch] = runkMeans1(X, initial_centroids, max_iters, false, kmeans_alpha);
        %[centroids_patch, idx_patch] = runkMeans(X, initial_centroids, max_iters, false);
        
        if (any(isnan(centroids_patch)) == true)
            idx_nan = isnan(centroids_patch);
            b = centroids_patch(~isnan(centroids_patch));
            centroids_patch(idx_nan) = mean(b);
        end
        
        context.centroids = centroids_patch;
        %centroids_patch_brightness = max(centroids_patch, [], 2);
        %centroids_patch_brightness = sum(centroids_patch.^2, 2).^.5;
        centroids_patch_brightness = exp(centroids_patch - additive);
%         centroids_patch_brightness = zeros(size(centroids_patch));
%         for k = 1:K
%             centroids_patch_brightness(k) = mean(X_hat(idx_patch==k));
%         end
        %centroids_patch_brightness = centroids_patch;

        C = sort(centroids_patch_brightness);
        I = zeros(size(C));
        for k = 1:K
            ii = find(centroids_patch_brightness(k) == C);
            try
                I(k) = ii;
            catch
                a = 55;
            end
        end
        centroids_adjustment = repmat((I - 0.5) / (K) ./ centroids_patch_brightness, 1, 3);
        centroids(n,:,:) = centroids_adjustment;
        idx(:,n) = idx_patch;
        
    end
end
% Image Compression
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
        max_adjustment_coeff = max((adjustment_coeff(:)));
        adjustment_coeff = adjustment_coeff/max_adjustment_coeff;
        %adjustment_coeff = bfilter2(adjustment_coeff,5,context.sigma);
        adjustment_coeff = imgaussfilt(adjustment_coeff, context.sigma);
        %adjustment_coeff = medfilt3(adjustment_coeff);
        adjustment_coeff = adjustment_coeff * max_adjustment_coeff;
        X_recovered = X_patch .* adjustment_coeff;
        X_recovered = reshape(X_recovered, image_patch(1), image_patch(2), 3);
        image_recovered(n,:,:,:) = X_recovered;
    end
end
% Putting image together
image_clustered = zeros(size(A));
n = 1;
for i = 1 : image_patch(1) : img_size(1) - image_patch(1) + 1
    for j = 1 : image_patch(2) : img_size(2) - image_patch(2) + 1
        image_clustered(i:i+image_patch(1)-1,j:j+image_patch(2)-1,:) = image_recovered(n, :, :, :);
        n = n + 1;
    end
end

out = image_clustered;

end

