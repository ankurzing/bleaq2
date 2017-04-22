function [eliteIndivLS, llEliteIndivLS,localSearch] = doLocalSearch2(archive, Indiv, ulDimMin, ulDimMax, llDimMin, llDimMax,flag,testProblemName)
    
    eliteIndiv = Indiv.upper; llEliteIndiv = Indiv.lower;
    ulDim = length(eliteIndiv); llDim = length(llEliteIndiv); 
    range = 0.05;
    
    doExactLSFlag = flag.localSearch;
    
    % obtain psi and phi mappings around eliteIndiv 
    [psiMapping,phiMapping] = getMappings(eliteIndiv,archive.tag1);
    functionLowerLevelVariables = psiMapping.function;
    functionLowerLevelObjective = phiMapping.function{:};
      
    localSearch.psiMSE = psiMapping.validMSE; localSearch.phiMSE = phiMapping.validMSE;
    
    if (psiMapping.validMSE <= phiMapping.validMSE) || flag.dontDoPhi == 1 || (rand <=0.5)
        if doExactLSFlag == 1
            localSearchMethod = 1; % psi w/ EXACT obj. func.
        else
            localSearchMethod = 2; % psi w/ Approx. obj. func.
        end
    else
        if doExactLSFlag == 1
            localSearchMethod = 3; % phi w/ EXACT obj. func.
        else
            localSearchMethod = 4; % phi w/ Approx. obj. func.
        end       
    end
    
    % data preparation
    archiveLS = [archive.tag1;archive.tag0];
    upper = cell2mat(cellfun(@(x) x.upper, archiveLS, 'UniformOutput',false));
    for j=1:size(upper,1)
        distances(j) = sum((eliteIndiv - upper(j,:)).^2);
    end
    [~, I] = sort(distances);
    
    if localSearchMethod == 2
        archiveConsidered = (ulDim+1)*(ulDim+2)/2+2*(ulDim)+ulDim;
        % if run psi w/ approx. obj. func. 
        % need to approximate F(x) and G(x) with Tag 1 member
        archivePsi.upper = upper(I(1:archiveConsidered),:);
        archivePsi.lower = cell2mat(cellfun(@(x) x.lower, archiveLS(I(1:archiveConsidered)), 'UniformOutput',false));
        archivePsi.functionValue = cell2mat(cellfun(@(x) x.functionValue, archiveLS(I(1:archiveConsidered)), 'UniformOutput',false));
        archivePsi.equalityConstrVals = cell2mat(cellfun(@(x) x.equalityConstrVals, archiveLS(I(1:archiveConsidered)), 'UniformOutput',false));
        archivePsi.inequalityConstrVals = cell2mat(cellfun(@(x) x.inequalityConstrVals, archiveLS(I(1:archiveConsidered)), 'UniformOutput',false));	

        approxPsi.function = quadApprox(archivePsi.functionValue, [archivePsi.upper]);
        if size(archivePsi.equalityConstrVals,2)~=0
	        for i=1:size(archivePsi.equalityConstrVals,2)  
                approxPsi.equalityConstr{i} = quadApprox(archivePsi.equalityConstrVals(:,i), [archivePsi.upper]);
	        end
	    else
	        approxPsi.equalityConstr = [];
        end
	    if size(archivePsi.inequalityConstrVals,2)~=0
	        for i=1:size(archivePsi.inequalityConstrVals,2)
                approxPsi.inequalityConstr{i} = quadApprox(archivePsi.inequalityConstrVals(:,i), [archivePsi.upper]);
	        end
	    else
	        approxPsi.inequalityConstr = [];
        end
    end
    
    if localSearchMethod == 4
        % if run phi w/ approx. obj. func. 
        % need to approximate F(x,y), G(x,y), f(x,y), g(x,y) w/ random
        % member
        archiveConsidered = (ulDim+llDim+1)*(ulDim+llDim+2)/2 + 2*(ulDim) + ulDim;
        
        archivePhi.upper = upper(I(1:archiveConsidered),:);
        archivePhi.lower = cell2mat(cellfun(@(x) x.lower, archiveLS(I(1:archiveConsidered)), 'UniformOutput',false));
        archivePhi.functionValue = cell2mat(cellfun(@(x) x.functionValue, archiveLS(I(1:archiveConsidered)), 'UniformOutput',false));
        archivePhi.equalityConstrVals = cell2mat(cellfun(@(x) x.equalityConstrVals, archiveLS(I(1:archiveConsidered)), 'UniformOutput',false));
        archivePhi.inequalityConstrVals = cell2mat(cellfun(@(x) x.inequalityConstrVals, archiveLS(I(1:archiveConsidered)), 'UniformOutput',false));	
        archivePhi.llFunctionValue = cell2mat(cellfun(@(x) x.llFunctionValue, archiveLS(I(1:archiveConsidered)), 'UniformOutput',false));
        archivePhi.llEqualityConstrVals = cell2mat(cellfun(@(x) x.llEqualityConstrVals, archiveLS(I(1:archiveConsidered)), 'UniformOutput',false));
        archivePhi.llInequalityConstrVals = cell2mat(cellfun(@(x) x.llInequalityConstrVals, archiveLS(I(1:archiveConsidered)), 'UniformOutput',false));
      
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

		% Phi-function with Approximated F(x,y) and f(x,y)
        % Needs to approximate f(x,y), g(x,y) & h(x,y)
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
    end
         

    if (localSearchMethod == 1)
        % Exact Psi-function approximation based local search    
        options = optimset('Algorithm','sqp','Display','off');
        [lb,ub] = createLocalSearchBound([eliteIndiv],[ulDimMin],[ulDimMax],range);
        [eliteIndivLS,~,EXITFLAG,OUTPUT] = fmincon(@(x) -approximatedFunctionPsi(x,llDimMin,llDimMax,functionLowerLevelVariables,testProblemName),[eliteIndiv],[],[],[],[],lb,ub,@(x) approximatedConstraintsPsi(x,llDimMin,llDimMax,functionLowerLevelVariables,testProblemName),options);
        llEliteIndivLS = llEliteIndiv;
        localSearch.method = 'Psi';
        localSearch.termination = EXITFLAG;
        localSearch.functionEvaluation = OUTPUT.funcCount;
        return;
    end

    if (localSearchMethod == 3)
    	% Exact Phi-function approximation based local search 
        options = optimset('Algorithm','sqp','Display','off');
        [lb,ub] = createLocalSearchBound([eliteIndiv llEliteIndiv],[ulDimMin llDimMin],[ulDimMax llDimMax],range);
        [eliteIndivFull,~,EXITFLAG,OUTPUT] = fmincon(@(x) -approximatedFunctionPhi(x,ulDim,testProblemName),[eliteIndiv llEliteIndiv],[],[],[],[],lb,ub,@(x) approximatedConstraintsPhi(x,functionLowerLevelObjective,ulDim, llDim,testProblemName),options);
        eliteIndivLS = eliteIndivFull(1:ulDim);
        llEliteIndivLS = eliteIndivFull(ulDim+1:end);
        localSearch.method = 'Phi';
        localSearch.termination = EXITFLAG;
        localSearch.functionEvaluation = OUTPUT.funcCount;
        return;
    end

    if (localSearchMethod == 2)
    	options = optimset('Algorithm','sqp','Display','off');      
        % psi-mapping based local search w/ approximated obj. func. 
        lb = ulDimMin; ub = ulDimMax;
        [eliteIndivLS,~,EXITFLAG,OUTPUT] = fmincon(@(x) -approximatedFunction(x,approxPsi.function),[eliteIndiv],[],[],[],[],lb,ub,@(x) approximatedConstraints(x,approxPsi.equalityConstr,approxPsi.inequalityConstr),options);
        llEliteIndivLS = llEliteIndiv;
        localSearch.method = 'Approx';
        localSearch.termination = EXITFLAG;
        localSearch.functionEvaluation = OUTPUT.funcCount;
        return;
    end
    
    if (localSearchMethod == 4)
        % phi-mapping based local search w/ approximated obj. func.
        options = optimset('Algorithm','sqp','Display','off');
        lb = [ulDimMin llDimMin]; ub = [ulDimMax llDimMax];
        [eliteIndivFull,~,EXITFLAG,OUTPUT] = fmincon(@(x) -approximatedFunction(x,approxPhi.function),[eliteIndiv llEliteIndiv],[],[],[],[],lb,ub,@(x) approximatedConstraintsPhi2(x,approxPhi.equalityConstr,approxPhi.inequalityConstr, approxPhi.llFunction, functionLowerLevelObjective, approxPhi.llEqualityConstr, approxPhi.llInequalityConstr, ulDim, llDim),options);
        eliteIndivLS = eliteIndivFull(1:ulDim);
        llEliteIndivLS = eliteIndivFull(ulDim+1:end);
        localSearch.method = 'ApproxPhi';
        localSearch.termination = EXITFLAG;
        localSearch.functionEvaluation = OUTPUT.funcCount;    
        return;
    end
return

function functionValue = approximatedFunctionPsi(xu,llDimMin,llDimMax,psiFunction,testProblemName)
    for j = 1:size(psiFunction,2)
        xl(j) = psiFunction{j}.constant + xu*psiFunction{j}.linear + xu*psiFunction{j}.sqmatrix*xu';
    end
    % check if predicted lower level optimal solution is outside the bound
    xl = checkLimits(xl,llDimMin,llDimMax);
	functionValue = ulTestProblem(xu, xl, testProblemName);

function functionValue = approximatedFunctionPhi(pop, ulDim, testProblemName)
    xu = pop(:,1:ulDim); xl = pop(:,ulDim+1:end);
	functionValue = ulTestProblem(xu,xl,testProblemName);
        
function approxFunctionValue = approximatedFunction(pop, parameters)
    approxFunctionValue = parameters.constant + pop*parameters.linear + pop*parameters.sqmatrix*pop';

function [c, ceq] = approximatedConstraintsPsi(xu,llDimMin,llDimMax,psiFunction, testProblemName)
    for j = 1:size(psiFunction,2)
        xl(j) = psiFunction{j}.constant + xu*psiFunction{j}.linear + xu*psiFunction{j}.sqmatrix*xu';
    end
    % check if predicted lower level optimal solution is outside the bound
    xl = checkLimits(xl,llDimMin,llDimMax);
    [~,ceq,c] = ulTestProblem(xu, xl, testProblemName);
    
function [c, ceq] = approximatedConstraintsPhi(pop, parametersPhiFunction, dimULPop, dimLLPop, testProblemName)

    ulPop = pop(:,1:dimULPop);
    llPop = pop(:,dimULPop+1:dimULPop+dimLLPop);
    
    [~,ceqUL,cUL] = ulTestProblem(ulPop, llPop, testProblemName);
    [c1,ceqLL,cLL] = llTestProblem(llPop, testProblemName, ulPop);
   
    ceq = [ceqUL ceqLL]; c = [cUL cLL];    
    
    n=length(c);
  
    c2 = (parametersPhiFunction.constant + ulPop*parametersPhiFunction.linear + ulPop*parametersPhiFunction.sqmatrix*ulPop');
    c(n+1) = c2 - c1;
    
function [c, ceq] = approximatedConstraints(pop, parametersEqualityConstr, parametersInequalityConstr)

    if ~isempty(parametersEqualityConstr)
        for i=1:length(parametersEqualityConstr)
            ceq(i) = parametersEqualityConstr{i}.constant + pop*parametersEqualityConstr{i}.linear + pop*parametersEqualityConstr{i}.sqmatrix*pop';
        end
    else
        ceq = [];
    end
    
    if ~isempty(parametersInequalityConstr)
        for i=1:length(parametersInequalityConstr)
            c(i) = parametersInequalityConstr{i}.constant + pop*parametersInequalityConstr{i}.linear + pop*parametersInequalityConstr{i}.sqmatrix*pop';
        end
    else
        c = []; %Negative value suggests that there is no active inequality
    end

function [c, ceq] = approximatedConstraintsPhi2(pop, parametersEqualityConstr, parametersInequalityConstr, parametersLowerLevelFunction, parametersPhiFunction, parametersLLEqualityConstr, parametersLLInequalityConstr, dimULPop, dimLLPop)

    ulPop = pop(:,1:dimULPop);
    llPop = pop(:,dimULPop+1:dimULPop+dimLLPop);
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

    
function offsprings=checkLimits(offsprings, DimMin, DimMax)
 
    numOffsprings = size(offsprings,1);
    dimMinMatrix = DimMin(ones(1,numOffsprings),:);
    offsprings(offsprings<dimMinMatrix)=dimMinMatrix(offsprings<dimMinMatrix);
    dimMaxMatrix = DimMax(ones(1,numOffsprings),:);
    offsprings(offsprings>dimMaxMatrix)=dimMaxMatrix(offsprings>dimMaxMatrix);
    
function ulPop=generateRandomSamples(ulPopSize,ulDim,ulDimMin,ulDimMax)

    for i=1:ulPopSize
        ulPop(i,:) = ulDimMin + rand(1, ulDim).*(ulDimMax-ulDimMin);
    end
    
function [LB,UB] = createLocalSearchBound(Indv,DimMin,DimMax,range)
	diff = range.*(DimMax - DimMin);
	LB = Indv - diff;
	LB = checkLimits(LB,DimMin,DimMax);
	UB = Indv + diff;
	UB = checkLimits(UB,DimMin,DimMax);
    