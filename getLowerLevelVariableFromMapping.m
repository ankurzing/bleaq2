function [offspringsLowerLevelVariables,optTypeOffsprings,sumMSE,validMSE] = getLowerLevelVariableFromMapping(offsprings,psiMapping,phiMapping,ulDim,llDim,archive)
    
    % minimum size required for phi-mapping based lower level optimal
    % variable retrieval
    minPhiSize = (ulDim+llDim+1)*(ulDim+llDim+2)/2 + 2*(ulDim) + ulDim;
    runPhi = 0;
    if (length(archive.tag1)+length(archive.tag0)) >= minPhiSize
        runPhi = ((psiMapping.sumMSE>=phiMapping.sumMSE) && (psiMapping.validMSE>=phiMapping.validMSE)) || (rand<0.25);
    end
    
	offspringsLowerLevelVariables = getLowerLevelVariableFromPsi(offsprings,psiMapping.function);
    optTypeOffsprings = 'Quadratic Approximation using Psi';
    sumMSE = psiMapping.sumMSE;
    validMSE = psiMapping.validMSE;
        
	if runPhi == 1
		offspringsLowerLevelVariables = getLowerLevelVariableFromPhi(offsprings,offspringsLowerLevelVariables,phiMapping.function{:},[archive.tag1;archive.tag0],minPhiSize);
        optTypeOffsprings = 'Quadratic Approximation using Phi';
        sumMSE = phiMapping.sumMSE;
        validMSE = phiMapping.validMSE;
    end
    
end	

function [indvLowerLevelVariables] = getLowerLevelVariableFromPsi(indv,psiFunction)

	llDim = length(psiFunction);
	indvLowerLevelVariables = zeros(1,llDim);
	for j = 1:llDim
		indvLowerLevelVariables(j) = psiFunction{j}.constant + indv*psiFunction{j}.linear + indv*psiFunction{j}.sqmatrix*indv';
	end	
end

function [indvLowerLevelVariables] = getLowerLevelVariableFromPhi(indv,indvLowerLevelVariables,phiFunction,archive,archiveSize)
	% make sure to use mixed archive (tag1 + tag0)
	% prepare data for phi-function based lower level variables optimization
	archivePhi.upper = cell2mat(cellfun(@(x) x.upper, archive, 'UniformOutput',false));

    for j=1:size(archivePhi.upper,1)
        distances(j) = sum((indv - archivePhi.upper(j,:)).^2);
    end
    [~, I] = sort(distances);
    I = I(1:archiveSize);
    archive = archive(I);

    archivePhi.upper = archivePhi.upper(I,:);
    archivePhi.lower = cell2mat(cellfun(@(x) x.lower, archive, 'UniformOutput',false));
	archivePhi.functionValue = cell2mat(cellfun(@(x) x.functionValue, archive, 'UniformOutput',false));
	archivePhi.equalityConstrVals = cell2mat(cellfun(@(x) x.equalityConstrVals, archive, 'UniformOutput',false));
	archivePhi.inequalityConstrVals = cell2mat(cellfun(@(x) x.inequalityConstrVals, archive, 'UniformOutput',false));	
	archivePhi.llFunctionValue = cell2mat(cellfun(@(x) x.llFunctionValue, archive, 'UniformOutput',false));
    archivePhi.llEqualityConstrVals = cell2mat(cellfun(@(x) x.llEqualityConstrVals, archive, 'UniformOutput',false));
    archivePhi.llInequalityConstrVals = cell2mat(cellfun(@(x) x.llInequalityConstrVals, archive, 'UniformOutput',false));
    
    approxPhi.function = quadApprox(archivePhi.functionValue, [archivePhi.upper archivePhi.lower]);
    if ~isempty(archivePhi.equalityConstrVals)
        for i=1:size(archivePhi.equalityConstrVals,2)		            
            approxPhi.equalityConstr{i} = quadApprox(archivePhi.equalityConstrVals(:,i), [archivePhi.upper archivePhi.lower]);
        end
    else
        approxPhi.equalityConstr = [];
    end
    if ~isempty(archivePhi.inequalityConstrVals)
        for i=1:size(archivePhi.inequalityConstrVals,2)
           approxPhi.inequalityConstr{i} = quadApprox(archivePhi.inequalityConstrVals(:,i), [archivePhi.upper archivePhi.lower]);            
        end
    else
        approxPhi.inequalityConstr = [];
    end

    approxPhi.llFunction = quadApprox(archivePhi.llFunctionValue, [archivePhi.upper archivePhi.lower]);

    if ~isempty(archivePhi.llEqualityConstrVals)
        for i=1:size(archivePhi.llEqualityConstrVals,2)
            approxPhi.llEqualityConstr{i} = quadApprox(archivePhi.llEqualityConstrVals(:,i), [archivePhi.upper archivePhi.lower]);
        end
    else
        approxPhi.llEqualityConstr = [];
    end

    if ~isempty(archivePhi.llInequalityConstrVals)
        for i=1:size(archivePhi.llInequalityConstrVals,2)
            approxPhi.llInequalityConstr{i} = quadApprox(archivePhi.llInequalityConstrVals(:,i), [archivePhi.upper archivePhi.lower]);
        end
    else
        approxPhi.llInequalityConstr = [];
    end
    lb = min(archivePhi.lower); ub = max(archivePhi.lower);
    options = optimset('Algorithm','sqp','Display','off');  
    [indvLowerLevelVariables,FVAL,EXITFLAG,OUTPUT] = fmincon(@(xl) -approximatedFunction(xl,indv,approxPhi.function),...
                    indvLowerLevelVariables,[],[],[],[],lb,ub,@(xl) approximatedConstraints(xl,indv,...
                    approxPhi.equalityConstr,approxPhi.inequalityConstr, approxPhi.llFunction,...
                    phiFunction, approxPhi.llEqualityConstr, approxPhi.llInequalityConstr),options);  
end

function approxFunctionValue = approximatedFunction(xl, xu, parameters)
    pop = [xu xl];
    approxFunctionValue = parameters.constant + pop*parameters.linear + pop*parameters.sqmatrix*pop'; 
end

function [c, ceq] = approximatedConstraints(llPop,ulPop, parametersEqualityConstr, parametersInequalityConstr, parametersLowerLevelFunction, parametersPhiFunction, parametersLLEqualityConstr, parametersLLInequalityConstr)

    pop = [ulPop,llPop];
    ceq = [];
    c = [];
    if ~isempty(parametersEqualityConstr)
        for i=1:length(parametersEqualityConstr)
            ceq(i) = parametersEqualityConstr{i}.constant + pop*parametersEqualityConstr{i}.linear + pop*parametersEqualityConstr{i}.sqmatrix*pop';
        end
    end
    n = length(ceq);
    if ~isempty(parametersLLEqualityConstr)
        for i=1:length(parametersLLEqualityConstr)
            ceq(n+i) = parametersLLEqualityConstr{i}.constant + pop*parametersLLEqualityConstr{i}.linear + pop*parametersLLEqualityConstr{i}.sqmatrix*pop';
        end
    end
    
    if ~isempty(parametersInequalityConstr)
        for i=1:length(parametersInequalityConstr)
            c(i) = parametersInequalityConstr{i}.constant + pop*parametersInequalityConstr{i}.linear + pop*parametersInequalityConstr{i}.sqmatrix*pop';
        end
    end
    n=length(c);
    if ~isempty(parametersLLInequalityConstr)
        for i=1:length(parametersLLInequalityConstr)
            c(i) = parametersLLInequalityConstr{i}.constant + pop*parametersLLInequalityConstr{i}.linear + pop*parametersLLInequalityConstr{i}.sqmatrix*pop';
        end
    end
    n=length(c);
    
    c1 = (parametersLowerLevelFunction.constant + pop*parametersLowerLevelFunction.linear + pop*parametersLowerLevelFunction.sqmatrix*pop');
    c2 = (parametersPhiFunction.constant + ulPop*parametersPhiFunction.linear + ulPop*parametersPhiFunction.sqmatrix*ulPop');
    c(n+1) = c2 - c1;
end
