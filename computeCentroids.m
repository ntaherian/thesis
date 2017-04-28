function centroids = computeCentroids(X, idx, K)


[m n] = size(X);

centroids = zeros(K, n);


for i = 1 : K
   points = find(idx == i);
   centroids(i,:) = (sum(X(points,:),1)) ./ size(points,1); 
end


end

