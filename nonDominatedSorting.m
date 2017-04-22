function [paretoRank] = nonDominatedSorting(functionValue, paretoRank)

noOfObj = size(functionValue,2);
popSize = size(functionValue,1);

if nargin==1
    paretoRank = zeros(popSize,1);
end

if popSize==0
    paretoRank = [];
    return
end

paretoRank = paretoRank+1;

for i=1:popSize
    dominance(i) = checkDominance(functionValue(i,:),functionValue);
end

[paretoRank(dominance)] = nonDominatedSorting(functionValue(dominance,:), paretoRank(dominance));


function dominance = checkDominance(objVectorToCheck, vectorsPool)
%Returns 1/true if objVectorToCheck is dominated by any member in
%vectorsPool else it returns 0/false
noOfObj = size(objVectorToCheck,2);
poolSize = size(vectorsPool,1);

for i=1:poolSize
    if objVectorToCheck~=vectorsPool(i,:)
        if objVectorToCheck<=vectorsPool(i,:)
            dominance = true;
            return;
        end
    end
end
dominance = false;
return;