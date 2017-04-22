function [crowdingDistance] = calculateCrowdingDistance(functionValue, paretoRank)

noOfObj = size(functionValue,2);
popSize = size(functionValue,1);

crowdingDistance = zeros(popSize,1);
for i=1:max(paretoRank)
    indices = find(paretoRank==i);
    [crowdingDistance(indices)] = calculateFrontCrowdingDistance(functionValue(indices,:));
end



function [crowdingDistance] = calculateFrontCrowdingDistance(functionValue)

noOfObj = size(functionValue,2);
popSize = size(functionValue,1);

sortingRank = zeros(size(functionValue));
sortedFunctionValue = zeros(size(functionValue));

for i=1:noOfObj
    [sortedFunctionValue(:,i), sortingIndices(:,i)] = sort(functionValue(:,i),'descend');
    [~,sortingRank(:,i)] = sort(sortingIndices(:,i),'ascend');
end

crowdingDistance = zeros(popSize,1);
for i=1:popSize
    for j=1:noOfObj
        lowerRank = find(sortingRank(i,j)+1==sortingRank(:,j));
        upperRank = find(sortingRank(i,j)-1==sortingRank(:,j));
        if ~isempty(lowerRank) && ~isempty(upperRank)
            crowdingDistance(i) = crowdingDistance(i) + abs(functionValue(i,j)-functionValue(lowerRank,j)) + abs(functionValue(i,j)-functionValue(upperRank,j));
        else
            crowdingDistance(i) = crowdingDistance(i)+Inf;
        end
    end
end