clc
clear all
close all

DataTable = readtable('house_prices_data_training_data.csv');
X = table2array(DataTable(1:17999,4:21));

% normalizing the X %
for i=1:length(X(1,:))
    X(:,i) = (X(:,i)-mean(X(:,i)))/std(X(:,i)); 
end

k = 10;
[ rows columns ] = size(X);
costVector = zeros(1,10);


for k = 1:8
    
    % initializing the centroids randomly %
    centroids = zeros( rows , columns);
    % initializing index %
    initial_index = randperm(rows);
    
    centroids = X(initial_index(1:k),:);
    prevCentroids = zeros(size(centroids));
    
    indices = zeros(size(X,1), 1);
    dist = zeros(rows,k);
    
    % flag to stop when prevcentroids is equal to centroids %
    flag = true;
    iter = 0;
    
    while(flag)
       for i = 1:rows
            for j = 1:k
                % get the distance between the data and centroid
                dist(i, j) = sum((X(i,:) - centroids(j, :)).^2);
            end
       end
        
        for i = 1:rows
            indices(i) = find(dist(i,:) == min(dist(i,:)));
        end
        
        for i = 1 : k
            
            clustering = X(find(indices == i), :);
            centroids(i, :) = mean(clustering);
            cost = 0;
            
             for z = 1 : size(clustering,1)
                cost = cost + sum((clustering(z,:) - centroids(i,:)).^2)/rows;
             end
            
            costVector(1,k) = cost;
            
        end
        
        if prevCentroids == centroids
            flag = false;
        end
        
        prevCentroids = centroids;
        
        iter = iter + 1;
        
    end
end
[ o bestKvalue ] = min(costVector);

noClusters = 1:10;

plot(noClusters, costVector)