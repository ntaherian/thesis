function idx = findClosestCentroids(X, centroids)

K = size(centroids, 1);

idx = zeros(size(X,1), 1);


distance_matrix = pdist2(X,centroids,'euclidean');
for i = 1:size(X,1)
    [~,index] = min(distance_matrix(i,:));
    idx(i,:) = index;
end

