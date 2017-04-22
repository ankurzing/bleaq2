function [psiMapping,phiMapping,lies] = getMappings(indv, archive, option)

    % construct psi-mapping or phi-mapping or both
    % indv - the point that local mappings are constructed around
    % archive - data archive structure (tag1 - optimal pair, tag0 - random pair)
    psiApprox = []; phiApprox = [];
    if nargin < 3
        option = 3;
    end
    constructPsi = 0; constructPhi = 0;
    if option == 1
        constructPsi = 1;
    elseif option == 2
        constructPhi = 1;
    elseif option == 3
        constructPsi = 1; constructPhi = 1;
    end
        
    ulDim = size(archive{1}.upper,2); llDim = size(archive{1}.lower,2); 
    
    %Minimum number of lower level points required for quadratic approximation
    %of both mappings
    numberLowerLevelPointsTrain = (ulDim+1)*(ulDim+2)/2+2*(ulDim);
    numberLowerLevelPointsEval = ulDim;
    minNumberLowerLevelPoints = numberLowerLevelPointsTrain + numberLowerLevelPointsEval;

    quadArchive.upper = cell2mat(cellfun(@(x) x.upper, archive, 'UniformOutput',false));
    quadArchive.lower = cell2mat(cellfun(@(x) x.lower, archive, 'UniformOutput',false));
    
    %Only the members close to the offspring are used for quadratic approximation
    distances = zeros(size(quadArchive.upper,1)); %Required correction: Line added on 28102014
    for k=1:length(quadArchive)
        distances(k) = sum((indv - quadArchive.upper(k,:)).^2);
    end
    [~, I] = sort(distances);
     %Note: minNumberLowerLevelPoints = numberLowerLevelPointsTrain + numberLowerLevelPointsEval
    I = I(1:minNumberLowerLevelPoints);           
    sizeI = length(I);
    permut = randperm(sizeI);
    J = permut(1:sizeI-numberLowerLevelPointsEval);
    quadApproxMembers = I(J);
    %I stores the members being considered from the archive
    %quadApproxMembers stores the members from I used for a quadratic approximation
    setDiff = setdiff(I,quadApproxMembers);
        
    if constructPsi == 1      
        psiApprox = cell(1,llDim); sumMSE = 0;
        for j = 1:llDim 
            psiApprox{j} = quadApprox(quadArchive.lower(quadApproxMembers,j), quadArchive.upper(quadApproxMembers,:));
            %Calculating sum of MSE for all lower level variable approximations
            sumMSE = sumMSE+psiApprox{j}.mseNorm;
        end
        predictedLowerLevelVariables = zeros(length(setDiff),llDim);
        for k = 1:length(setDiff)
            for j=1:llDim
                predictedLowerLevelVariables(k,j) = psiApprox{j}.constant + quadArchive.upper(setDiff(k),:)*psiApprox{j}.linear + quadArchive.upper(setDiff(k),:)*psiApprox{j}.sqmatrix*quadArchive.upper(setDiff(k),:)';
            end
        end
        psiMapping.function = psiApprox;
        psiMapping.sumMSE = sumMSE;
        psiMapping.validMSE = mean(mean((predictedLowerLevelVariables - quadArchive.lower(setDiff,:)).^2));
    end
    
    if constructPhi == 1 
        quadArchive.llFunctionValue = cell2mat(cellfun(@(x) x.llFunctionValue, archive,'UniformOutput',false));
        llFunctionDim = size(quadArchive.llFunctionValue,2);
        sumMSE = 0; phiApprox = cell(1,llFunctionDim);
        for j = 1:llFunctionDim
            phiApprox{j} = quadApprox(quadArchive.llFunctionValue(quadApproxMembers,j), quadArchive.upper(quadApproxMembers,:));
            sumMSE = sumMSE+phiApprox{j}.mseNorm;
        end
        predictedLowerLevelObjective = zeros(length(setDiff),llFunctionDim);  
        for k = 1:length(setDiff)
            for j = 1:llFunctionDim
                predictedLowerLevelObjective(k,j) = phiApprox{j}.constant + quadArchive.upper(setDiff(k),:)*phiApprox{j}.linear + quadArchive.upper(setDiff(k),:)*phiApprox{j}.sqmatrix*quadArchive.upper(setDiff(k),:)';
            end
        end
        phiMapping.function = phiApprox;
        phiMapping.sumMSE = sumMSE;
        phiMapping.validMSE = mean(mean((predictedLowerLevelObjective - quadArchive.llFunctionValue(setDiff,:)).^2));
    end
    
    %Checks if the offspring lies in between the training points
    %and not outside
    lies = 1;
    for j=1:size(quadArchive.upper,2)
        if indv(j)>=max(quadArchive.upper(quadApproxMembers,j)) || indv(j)<=min(quadArchive.upper(quadApproxMembers,j))
            lies=0;
        end
    end
end