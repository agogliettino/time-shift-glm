function [distToNearNeighb,nullDist,pval] = findNearestNeighbor_shuffleLabels(var1,var2,labels,numNeigb,numIter)
% labels is a vector of 0 and 1's - within group and not within group
% var1 is on the x-axis, var2 is on the y-axis
% numIter is the number of shuffles


labels = logical(labels);

% make sure these are column vectors
var1 = reshape(var1,numel(var1),1);
var2 = reshape(var2,numel(var2),1);

% find the nearest neigbhor for actual labels
var1_label = var1(labels);
var2_label = var2(labels);
varMat = [var1_label var2_label];

% compute mean distance to k nearest neighbors
distance = pdist2(varMat,varMat,'euclidean','Smallest',numNeigb+1);
distToNearNeighb = mean(mean(distance(2:end,:)));

nullDist = nan(numIter,1);

if isfinite(distToNearNeighb)
    
    for k = 1:numIter
        labels_new = labels(randperm(numel(labels)));
        var1_label = var1(labels_new);
        var2_label = var2(labels_new);
        varMat = [var1_label var2_label];
        distance = pdist2(varMat,varMat,'euclidean','Smallest',numNeigb+1);
        nullDist(k) = mean(mean(distance(2:end,:)));
    end
end

pval = sum(distToNearNeighb < nullDist)/numIter;


return