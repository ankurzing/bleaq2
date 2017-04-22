function [eliteIndiv, eliteFunctionValue, eliteConstr, tag, freeVariablesYesNo, gen]=llSearch(testProblemName, llPopSize, llMaxGens, llDim, llDimMin, llDimMax,  upperLevelVariables, llEpsilonStopping, llMember, llMemberVariance, llPop)

global testProblemNameLowerLevel;
global ulMember;

testProblemNameLowerLevel=testProblemName;
ulMember=upperLevelVariables;


%Call as [eliteFitness, eliteIndiv]=g3pcx(testProblemNameLowerLevel, popSize, maxGens, dim, dimMin, dimMax)
%Algorithm for maximization problem
%feasibility g(x)<0 and h(x)=0

if ~isempty(llPop) 
    popSize=size(llPop,1);         % The size of the population
else
    popSize=llPopSize;         % The size of the population
end
maxGens=llMaxGens;         % The maximum number of generations allowed in a run
dim=llDim;                 % Number of dimensions
dimMin=llDimMin;           % Minimum value accross dimensions
dimMax=llDimMax;           % Maximum value accross dimensions

probCrossover=0.9;         % The probability of crossover. 
probMutation=1/llDim;      % The mutation probability
%probMutation=0.1;          % The mutation probability

crossoverType=1;           % 0 => no crossover
                           % 1 => planar pcx crossover
                           % 2 => normal pcx crossover
                           % 3 => ucx crossover
                           % 4 => adaptive planar pcx crossover

mutationType=1;            % 0 => no mutation
                           % 1 => polynomial mutation
                           % 2 => variance based polynomial mutation

visualizationFlag=0;       % 0 => don't visualize
                           % 1 => visualize

verboseFlag=0;             % 0 => run quietly
                           % 1 => display details of each generation

numParents=3;
numOffsprings=2;
kappa = popSize;

violationWindowEquality = 1e-1;
violationWindowInequality = 1e-6;
violationWindowConstraint = 1e-6;

linearApproxError = 1e-14;
epsilonZero = 1e-6;
%llEpsilonStopping = 1e-5;
%llEpsilonStoppingImprovement = llEpsilonStopping;
%llImprovementGenDiff = 500;
llEpsilonStoppingImprovement = 1e-4;
improvementGenDiff = 100;
times = 1;
certainlyDoGaSearch = 0;

eliteIndiv=[];
eliteFitness=-realmax;

%Initialization
allPop = [];
allFunctionValue = [];
allEqualityConstrVals = [];
allInequalityConstrVals = [];
tag = 0;
gen = 0;

for iter=1:times
    popSize = (dim+1)*(dim+2)/2+dim;
    if ~isempty(llMember)
        for i=1:popSize-1
            if ~isempty(llMemberVariance)
                %Following line added on 21-08-2013
                llMemberVariance(llMemberVariance<epsilonZero)=epsilonZero;
                pop_new(i,:)=llMember + (epsilonZero+sqrt(llMemberVariance)).*randn(size(llMember));
            else
                pop_new(i,:)=llMember + 0.1/times*(dimMax-dimMin).*randn(size(llMember));
            end
        end
        pop_new(popSize,:) = llMember;
        pop_new=checkLimitsReflection(pop_new, dimMin, dimMax);
    else
        pop_new=initialize(popSize,dim,dimMin,dimMax);
        llMember = pop_new(popSize,:);
    end

    [functionValue_new equalityConstrVals_new inequalityConstrVals_new]=llTestProblem(pop_new, testProblemName, ulMember);

    if isempty(allPop)
        allPop = pop_new;
        allFunctionValue = functionValue_new;
        allEqualityConstrVals = equalityConstrVals_new;
        allInequalityConstrVals = inequalityConstrVals_new;
    else
        allPop(end+1:end+size(pop_new,1),:) = pop_new;
        allFunctionValue(end+1:end+size(functionValue_new,1)) = functionValue_new;
        if ~isempty(equalityConstrVals_new)
            allEqualityConstrVals(end+1:end+size(equalityConstrVals_new,1),:) = equalityConstrVals_new;
        else
            allEqualityConstrVals = [];
        end
        if ~isempty(inequalityConstrVals_new)
            allInequalityConstrVals(end+1:end+size(inequalityConstrVals_new,1),:) = inequalityConstrVals_new;
        else
            allInequalityConstrVals = [];
        end
    end

    pop = allPop;
    functionValue = allFunctionValue;
    equalityConstrVals = allEqualityConstrVals;
    inequalityConstrVals = allInequalityConstrVals;

    numEqualityConstr = size(equalityConstrVals,2);
    numInequalityConstr = size(inequalityConstrVals,2);
    
	approx.function = quadApprox(functionValue, pop);
        
    if size(equalityConstrVals,2)~=0
        for i=1:size(equalityConstrVals,2)
            approx.equalityConstr{i} = quadApprox(equalityConstrVals(:,i), pop);
        end
    else
        approx.equalityConstr = [];
    end

    if size(inequalityConstrVals,2)~=0
        for i=1:size(inequalityConstrVals,2)
            approx.inequalityConstr{i} = quadApprox(inequalityConstrVals(:,i), pop);
        end
    else
        approx.inequalityConstr = [];
    end

    options = optimset('Algorithm','active-set'); % run active-set algorithm
    options = optimset('Display','off','TolX',1e-10,'TolFun',1e-10);
    
    if isLinear(approx,epsilonZero)
        if ~isempty(approx.equalityConstr)
            for i=1:length(approx.equalityConstr)
                A_equality(i,:) = approx.equalityConstr{i}.linear;
                b_equality(i) = -approx.equalityConstr{i}.constant;
            end
        else
            A_equality = [];
            b_equality = [];
        end
        if ~isempty(approx.inequalityConstr)
            for i=1:length(approx.inequalityConstr)
                A_inequality(i,:) = approx.inequalityConstr{i}.linear;
                b_inequality(i) = -approx.inequalityConstr{i}.constant;
            end
        else
            A_inequality = [];
            b_inequality = [];
        end
        
        optionsLinprog = optimset('Display','off');
        [eliteIndiv,FVAL,EXITFLAG,OUTPUT,LAMBDA] = linprog(-approx.function.linear,A_inequality,b_inequality',A_equality,b_equality',dimMin,dimMax,llMember,optionsLinprog);
        eliteIndiv = eliteIndiv';
        LAMBDA.eqnonlin = [];
        LAMBDA.ineqnonlin = [];
    else
        [eliteIndiv,FVAL,EXITFLAG,OUTPUT,LAMBDA] = fmincon(@(x) -approximatedFunction(x,approx.function),llMember,[],[],[],[],dimMin,dimMax,@(x) approximatedConstraints(x,approx.equalityConstr,approx.inequalityConstr),options);
    end

    %Rows corresponding to bound constraints are not required in bordered Hessian matrix

    if size(equalityConstrVals,2)~=0
        for i=1:size(equalityConstrVals,2)
            gradH(i,:) =  eliteIndiv * 2*eye(dim,dim)*approx.equalityConstr{i}.sqmatrix + approx.equalityConstr{i}.linear';
        end
    else
        gradH = [];
    end

    if size(inequalityConstrVals,2)~=0
        for i=1:size(inequalityConstrVals,2)
            gradG(i,:) = eliteIndiv * 2*eye(dim,dim)*approx.inequalityConstr{i}.sqmatrix + approx.inequalityConstr{i}.linear';
        end
    else
        gradG = [];
    end

    %Note: Border rows for upper and lower bound constraints are not required

    bordered_hessian = zeros(numEqualityConstr+numInequalityConstr+dim,numEqualityConstr+numInequalityConstr+dim);

    if size(equalityConstrVals,2)~=0
        tempGradH = gradH;
        tempGradH(abs(LAMBDA.eqnonlin)<epsilonZero,:) = 0;
        bordered_hessian(1:numEqualityConstr,numEqualityConstr+numInequalityConstr+1:end) = tempGradH;
        bordered_hessian(numEqualityConstr+numInequalityConstr+1:end,1:numEqualityConstr) = tempGradH';
    end
    if size(inequalityConstrVals,2)~=0
        tempGradG = gradG;
        tempGradG(abs(LAMBDA.ineqnonlin)<epsilonZero,:) = 0;
        bordered_hessian(1+numEqualityConstr:numEqualityConstr+numInequalityConstr,numEqualityConstr+numInequalityConstr+1:end) = tempGradG;
        bordered_hessian(numEqualityConstr+numInequalityConstr+1:end,1+numEqualityConstr:numEqualityConstr+numInequalityConstr) = tempGradG';
    end
    bordered_hessian(end-dim+1:end,end-dim+1:end) = approx.function.sqmatrix;

    [R,jb] = rref(bordered_hessian,epsilonZero);
    boundVariables = jb - (numEqualityConstr+numInequalityConstr);

    k = 0;
    freeVariables = [];
    freeVariablesYesNo = zeros(1,llDim);
    %Check: "abs(eliteIndiv(i)-dimMin(i))>epsilonZero & abs(dimMax(i)-eliteIndiv(i))>epsilonZero" to ensure that variables on bounds are not not considered as free variables
    for i=1:dim
        if isempty(find(boundVariables==i)) & abs(eliteIndiv(i)-dimMin(i))>epsilonZero & abs(dimMax(i)-eliteIndiv(i))>epsilonZero
            k = k + 1;
            freeVariables(k) = i;
        end
    end

    if ~isempty(freeVariables)
        if mutationType == 1
            eliteIndiv(freeVariables) = mutation(llMember(freeVariables), probMutation);
        elseif mutationType == 2
            if ~isempty(llMemberVariance)
                eliteIndiv(freeVariables) = mutationVarianceBased(llMember(freeVariables), probMutation, llMemberVariance(freeVariables));
            else
                eliteIndiv(freeVariables) = mutationVarianceBased(llMember(freeVariables), probMutation, []);
            end
        end
        for i=1:size(freeVariables,2)
            j=size(approx.equalityConstr,2)+1;
            approx.equalityConstr{j}.constant = -eliteIndiv(freeVariables(i));
            approx.equalityConstr{j}.linear = zeros(dim,1);
            approx.equalityConstr{j}.linear(freeVariables(i)) = 1;
            approx.equalityConstr{j}.sqmatrix = 0;
        end
        options = optimset('Display','off','TolX',1e-10,'TolFun',1e-10);
        [eliteIndiv,FVAL,EXITFLAG,OUTPUT,LAMBDA] = fmincon(@(x) -approximatedFunction(x,approx.function),eliteIndiv,[],[],[],[],dimMin,dimMax,@(x) approximatedConstraints(x,approx.equalityConstr,approx.inequalityConstr),options);
        for i=1:length(freeVariables)
           freeVariablesYesNo(freeVariables(i)) = 1;
        end
        approx.equalityConstr = approx.equalityConstr(1:size(equalityConstrVals,2));
    end

    isBoundFeasible = false;
    if eliteIndiv>dimMin - epsilonZero
        if eliteIndiv<dimMax + epsilonZero
            isBoundFeasible = true;
        end
    end
    eliteIndiv(eliteIndiv<dimMin) = dimMin(eliteIndiv<dimMin);
    eliteIndiv(eliteIndiv>dimMax) = dimMax(eliteIndiv>dimMax);
    [eliteFunctionValue, eliteEqualityConstrVals, eliteInequalityConstrVals]=llTestProblem(eliteIndiv, testProblemName, ulMember);

    allPop(end+1,:) = eliteIndiv;
    allFunctionValue(end+1) = eliteFunctionValue;
    if ~isempty(eliteEqualityConstrVals)
        allEqualityConstrVals(end+1,:) = eliteEqualityConstrVals;
    else
        allEqualityConstrVals = [];
    end
    if ~isempty(eliteInequalityConstrVals)
        allInequalityConstrVals(end+1,:) = eliteInequalityConstrVals;
    else
        allInequalityConstrVals = [];
    end
    [allFitness, allConstrViolation] = assignFitness(allFunctionValue, allEqualityConstrVals, allInequalityConstrVals, violationWindowEquality, violationWindowInequality);
    if ~isempty(allConstrViolation)
        eliteConstrViolation = allConstrViolation(end);
    else
        eliteConstrViolation = [];
    end
    
    f=approximatedFunction(eliteIndiv,approx.function);
    [c, ceq]=approximatedConstraints(eliteIndiv,approx.equalityConstr,approx.inequalityConstr);
    d = sqrt((f-eliteFunctionValue)^2+sum((c-eliteInequalityConstrVals).^2)+sum((ceq-eliteEqualityConstrVals).^2));
    if d>epsilonZero || ~isFeasible(eliteIndiv,eliteConstrViolation,dimMin,dimMax,violationWindowConstraint) || ~isBoundFeasible
        llMember = eliteIndiv;
        if d<epsilonZero
            display('No lower level feasible solution found for given x_u.')          
            tag = 0;
            eliteConstr.constrViolation = eliteConstrViolation; 
            eliteConstr.equalityConstrVals = eliteEqualityConstrVals;
            eliteConstr.inequalityConstrVals = eliteInequalityConstrVals;
            return;
        end
        display('SQP unsuccessful at lower level.')
    else
        if certainlyDoGaSearch == 1
            tag = 0;
        else
            tag = 1;
        end
        break;
    end
end
[allFitness, allConstrViolation] = assignFitness(allFunctionValue, allEqualityConstrVals, allInequalityConstrVals, violationWindowEquality, violationWindowInequality);
[maxFitness, I] = max(allFitness);
eliteIndiv = allPop(I,:);
eliteFunctionValue = allFunctionValue(I);
if ~isempty(allConstrViolation)
    eliteConstrViolation = allConstrViolation(I);
else
    eliteConstrViolation = [];
end

% % added on 2017/1/22
% tag = 1;
if tag==1 
    if isnan(eliteFunctionValue)
        disp('Function value is NaN in lower level quadratic programming.');
        clear eliteIndiv;
        return;
    end
    display('SQP successful at lower level.');
    eliteConstr.constrViolation = eliteConstrViolation; 
    eliteConstr.equalityConstrVals = eliteEqualityConstrVals;
    eliteConstr.inequalityConstrVals = eliteInequalityConstrVals;
    return;
end
eliteIndivQuadratic = eliteIndiv;
eliteFunctionValueQuadratic = eliteFunctionValue;
eliteEqualityConstrValsQuadratic = eliteEqualityConstrVals;
eliteInequalityConstrValsQuadratic = eliteInequalityConstrVals;


if tag==0
    display('Doing lower level GA search.');
    if ~isempty(llPop) 
        if size(llPop,1)>llPopSize
            popSize=llPopSize;     % The size of the population
            r = randperm(size(llPop,1));
            r = r(1:llPopSize);
            pop=llPop(r,:);
        else
            popSize=llPopSize;
            pop=initialize(popSize,dim,dimMin,dimMax);
            %Replace some of the population members with members that have been
            %passed from the upper level
            pop(1:size(llPop,1),:)=llPop;
        end
    else
        popSize=llPopSize;         % The size of the population
        pop=initialize(popSize,dim,dimMin,dimMax);
    end
end
% evaluate the fitness of the population. The vector of fitness values 
% returned  must be of dimensions 1 x popSize.
tag = 1;
gen = 0;
[functionValue equalityConstrVals inequalityConstrVals]=llTestProblem(pop, testProblemName, ulMember);
if exist('eliteIndivQuadratic','var')
    popSize = popSize + 1;
    pop(popSize,:)=eliteIndivQuadratic;
    functionValue(popSize) = eliteFunctionValueQuadratic;
    if ~isempty(equalityConstrVals)
        equalityConstrVals(popSize, :) = eliteEqualityConstrValsQuadratic;
    end
    if ~isempty(inequalityConstrVals)
        inequalityConstrVals(popSize, :) = eliteInequalityConstrValsQuadratic;
    end
end
[fitnessVals constrViolation] = assignFitness(functionValue, equalityConstrVals, inequalityConstrVals, violationWindowEquality*(1+exp(-gen*numOffsprings/popSize)), violationWindowInequality);

alpha_not = sum(var(pop));
alpha = 1;

for gen=1:maxGens

    param.alpha = alpha;
    param.kappa = kappa;
    param.functionValue = functionValue;
    param.constrViolation = constrViolation;
    param.popSize = popSize;
    
%     if ~isempty(llPop)
%         [indicesParents fitnessParents functionValueParents constrViolationParents equalityConstrValsParents inequalityConstrValsParents parents]=selectParents2(pop, fitnessVals, functionValue, constrViolation, equalityConstrVals, inequalityConstrVals, numParents);
%     else
%         [indicesParents fitnessParents functionValueParents constrViolationParents equalityConstrValsParents inequalityConstrValsParents parents]=tournamentSelectSortedParents(pop, fitnessVals, functionValue, constrViolation, equalityConstrVals, inequalityConstrVals, numParents);
%     end
    [indicesParents fitnessParents functionValueParents constrViolationParents equalityConstrValsParents inequalityConstrValsParents parents]=selectParents2(pop, fitnessVals, functionValue, constrViolation, equalityConstrVals, inequalityConstrVals, numParents);
        
    
    if crossoverType == 1 
        %planar pcx crossover
        offsprings=planar_pcx_crossover(parents, probCrossover, numOffsprings, param);
        offsprings=checkLimits(offsprings, dimMin, dimMax);
    elseif crossoverType == 2 
        %normal pcx crossover
        offsprings=pcx_crossover(parents, probCrossover, numOffsprings, param);
        offsprings=checkLimits(offsprings, dimMin, dimMax);
    elseif crossoverType == 3
        %Uniform centric crossover
        offsprings=ucx_crossover(parents, probCrossover, numOffsprings, param);
        offsprings=checkLimits(offsprings, dimMin, dimMax);
    elseif crossoverType == 4
        %Adaptive planar pcx crossover
        offsprings=adaptive_planar_pcx_crossover(parents, probCrossover, numOffsprings, param);
        offsprings=checkLimits(offsprings, dimMin, dimMax);
    end
    
    if mutationType == 1 
        %polynomial mutation
        offsprings=mutation(offsprings, probMutation);
        offsprings=checkLimits(offsprings, dimMin, dimMax);
    elseif mutationType == 2
        %variance based polynomial mutation
        offsprings=mutationVarianceBased(offsprings, probMutation, var(pop));
        offsprings=checkLimits(offsprings, dimMin, dimMax);
    end
    
    [functionValueOffsprings, equalityConstrValsOffsprings, inequalityConstrValsOffsprings]=llTestProblem(offsprings, testProblemName, ulMember);
    if sum(isnan(functionValueOffsprings))>0
        disp('Function value is NaN in GA at lower level.');
        clear eliteIndiv;
        return;
    end
    [tempFitnessVals tempConstrViolation] = assignFitness([functionValue; functionValueOffsprings], [equalityConstrVals; equalityConstrValsOffsprings], [inequalityConstrVals; inequalityConstrValsOffsprings], violationWindowEquality*(1+exp(-gen*numOffsprings/popSize)), violationWindowInequality);
    fitnessVals = tempFitnessVals(1:popSize);
    fitnessOffsprings = tempFitnessVals(1+popSize:end);
    if ~isempty(tempConstrViolation)
        constrViolation = tempConstrViolation(1:popSize);
        constrViolationOffsprings = tempConstrViolation(1+popSize:end);
    else
        constrViolation = [];
        constrViolationOffsprings = [];
    end
    fitnessParents = fitnessVals(indicesParents);
    [pop fitnessVals functionValue constrViolation equalityConstrVals inequalityConstrVals] = update(pop, fitnessVals, functionValue, constrViolation, equalityConstrVals, inequalityConstrVals, indicesParents, parents, fitnessParents, functionValueParents, constrViolationParents, equalityConstrValsParents, inequalityConstrValsParents, offsprings, fitnessOffsprings, functionValueOffsprings, constrViolationOffsprings,  equalityConstrValsOffsprings, inequalityConstrValsOffsprings);
    
    [eliteFitness, maxIndex]=max(fitnessVals);
    eliteIndiv=pop(maxIndex,:);
    if ~isempty(constrViolation)
        eliteConstrViolation=constrViolation(maxIndex);
    else
        eliteConstrViolation=0;
    end
    eliteFunctionValue=functionValue(maxIndex);
    if ~isempty(equalityConstrVals)
        eliteEqualityConstrVals = equalityConstrVals(maxIndex,:);
    end
    if ~isempty(inequalityConstrVals)
        eliteInequalityConstrVals = inequalityConstrVals(maxIndex,:);
    end
    % added on 2016/08/16 to output equality & inequality constraint values
    % not just violation, which are used in doLocalSearch2 in ulSearch.m
    eliteConstr.constrViolation = eliteConstrViolation; 
    eliteConstr.equalityConstrVals = eliteEqualityConstrVals;
    eliteConstr.inequalityConstrVals = eliteInequalityConstrVals;
    % end 
    
    avgFitness = mean(fitnessVals);
    alpha = sum(var(pop))/alpha_not;
    if alpha>1
        alpha = 1;
    end
    
    eliteFunctionValueAtGen(gen) = eliteFunctionValue;
    if (gen>improvementGenDiff)
        if abs(eliteFunctionValueAtGen(gen)-eliteFunctionValueAtGen(gen-improvementGenDiff)) == 0
            beta = 0;
        elseif (abs(eliteFunctionValueAtGen(gen))+abs(eliteFunctionValueAtGen(1))) == 0
            beta = abs(eliteFunctionValueAtGen(gen)-eliteFunctionValueAtGen(gen-improvementGenDiff));
        else
            beta = abs(eliteFunctionValueAtGen(gen)-eliteFunctionValueAtGen(gen-improvementGenDiff))/(abs(eliteFunctionValueAtGen(gen))+abs(eliteFunctionValueAtGen(1)));
        end
    else
        beta = Inf;
    end
    
    % display the generation number, the average Fitness of the population,
    % and the maximum function value in the population
    if verboseFlag
        if constrViolation > 0
            display(['gen=' num2str(gen,'%.3d') '   At lower level no feasible solution found.']);
        else
            display(['At lower level gen=' num2str(gen,'%.3d') '   avgFitness=' ...
            num2str(avgFitness,'%3.3f') '   maxFitnessValue=' ...
            num2str(eliteFitness,'%3.3f')]);
        end
    end
    % Perform visualization
    if visualizationFlag
        figure(1)
        set (gcf, 'color', 'w');
        hold off
        plot(1:popSize,fitnessVals, '.');
        title(['Generation = ' num2str(gen) ', Average Fitness = ' sprintf('%0.3f', avgFitness)]);
        ylabel('Fitness');
        xlabel('Population members');
        drawnow;
    end
    
    if alpha<llEpsilonStopping & eliteConstrViolation<=0
        display('Lower level search done. Variance based termination condition met.');
%         [eliteIndivLS, eliteFunctionValueLS, eliteConstrViolationLS, tag]=lowerLevelLocalSearch(testProblemName, llDim, llDimMin, llDimMax,  upperLevelVariables, pop, functionValue, inequalityConstrVals, equalityConstrVals, eliteIndiv);
%         if tag==1
%             eliteIndivLS = eliteIndiv;
%             eliteFunctionValueLS = eliteFunctionValue;
%             eliteConstrViolationLS = eliteConstrViolation;
%         end
        return;
    end
    
    if beta<llEpsilonStoppingImprovement & eliteConstrViolation<=0
        display('Lower level search done. Improvement based termination condition met.');
%         [eliteIndivLS, eliteFunctionValueLS, eliteConstrViolationLS, tag]=lowerLevelLocalSearch(testProblemName, llDim, llDimMin, llDimMax,  upperLevelVariables, pop, functionValue, inequalityConstrVals, equalityConstrVals, eliteIndiv);
%         if tag==1
%             eliteIndivLS = eliteIndiv;
%             eliteFunctionValueLS = eliteFunctionValue;
%             eliteConstrViolationLS = eliteConstrViolation;
%         end
        return;
    end
    
end
tag = 0;
if eliteConstrViolation>0
    display('No feasible solution found at lower level.')
end
display('Lower level search done. Termination condition not met. Lower level solution might be inaccurate. Consider increasing lower level generations.');

function pop=initialize(popSize,dim,dimMin,dimMax)

    for i=1:popSize
        pop(i,:) = dimMin + rand(1, dim).*(dimMax-dimMin);
    end

function [indicesParents fitnessParents functionValueParents constrViolationParents equalityConstrValsParents inequalityConstrValsParents parents]=tournamentSelectParents(pop, fitnessVals, functionValue, constrViolation, equalityConstrVals, inequalityConstrVals, numParents)

    popSize = size(pop, 1);
    
    if ~isempty(constrViolation)
        t = find(constrViolation==0);
        if size(t,1)>2*numParents
            permut = randperm(size(t,1));
            s = t(permut(1:2*numParents));
        else
            [sortedFitness, I] = sort(fitnessVals, 'descend');
            s = I(1:2*numParents);
        end
    else
        permut = randperm(popSize);
        s = permut(1:2*numParents);
    end
    
    for i=1:2:2*numParents
        if fitnessVals(s(i))>fitnessVals(s(i+1))
            r((i+1)/2) = s(i);
        else
            r((i+1)/2) = s(i+1);
        end
    end
        
    parents = [pop(r,:)];
    fitnessParents = [fitnessVals(r)];
    functionValueParents = [functionValue(r)];

    if ~isempty(equalityConstrVals)
        equalityConstrValsParents = [equalityConstrVals(r,:)];
    else
        equalityConstrValsParents = [];
    end
    if ~isempty(inequalityConstrVals)
        inequalityConstrValsParents = [inequalityConstrVals(r,:)];
    else
        inequalityConstrValsParents = [];
    end
    if ~isempty(constrViolation)
        constrViolationParents = [constrViolation(r,:)];
    else
        constrViolationParents = [];
    end
    indicesParents = r;

function [indicesParents fitnessParents functionValueParents constrViolationParents equalityConstrValsParents inequalityConstrValsParents parents]=tournamentSelectSortedParents(pop, fitnessVals, functionValue, constrViolation, equalityConstrVals, inequalityConstrVals, numParents)

    %In this tournament selection, tournament is performed between two
    %random members and one is chosen as parent. The winners of the 
    %tournament are sorted by fitness and then returned.
    popSize = size(pop, 1);
    
    if ~isempty(constrViolation)
        t = find(constrViolation==0);
        if size(t,1)>2*numParents
            permut = randperm(size(t,1));
            s = t(permut(1:2*numParents));
        else
            [sortedFitness, I] = sort(fitnessVals, 'descend');
            s = I(1:2*numParents);
        end
    else
        permut = randperm(popSize);
        s = permut(1:2*numParents);
    end
    
    for i=1:2:2*numParents
        if fitnessVals(s(i))>fitnessVals(s(i+1))
            r((i+1)/2) = s(i);
        else
            r((i+1)/2) = s(i+1);
        end
    end
    
    [sortedFitnessParents, I] = sort(fitnessVals(r),'descend');
    parents = [pop(r(I),:)];
    fitnessParents = [fitnessVals(r(I))];
    functionValueParents = [functionValue(r(I))];

    if ~isempty(equalityConstrVals)
        equalityConstrValsParents = [equalityConstrVals(r(I),:)];
    else
        equalityConstrValsParents = [];
    end
    if ~isempty(inequalityConstrVals)
        inequalityConstrValsParents = [inequalityConstrVals(r(I),:)];
    else
        inequalityConstrValsParents = [];
    end
    if ~isempty(constrViolation)
        constrViolationParents = [constrViolation(r(I),:)];
    else
        constrViolationParents = [];
    end
    indicesParents = r(I);
    
function [indicesParents fitnessParents functionValueParents constrViolationParents equalityConstrValsParents inequalityConstrValsParents parents]=selectParents(pop, fitnessVals, functionValue, constrViolation, equalityConstrVals, inequalityConstrVals, numParents)

    popSize = size(pop, 1);
    r = randint(1,numParents-1,[1 popSize]);
    [eliteFitness, maxIndex]=max(fitnessVals);
    eliteIndiv=pop(maxIndex,:);
    
    parents = [pop(maxIndex,:); pop(r,:)];
    fitnessParents = [fitnessVals(maxIndex); fitnessVals(r)];
    functionValueParents = [functionValue(maxIndex); functionValue(r)];
    if ~isempty(equalityConstrVals)
        equalityConstrValsParents = [equalityConstrVals(maxIndex,:); equalityConstrVals(r,:)];
    else
        equalityConstrValsParents = [];
    end
    if ~isempty(inequalityConstrVals)
        inequalityConstrValsParents = [inequalityConstrVals(maxIndex,:); inequalityConstrVals(r,:)];
    else
        inequalityConstrValsParents = [];
    end
    if ~isempty(constrViolation)
        constrViolationParents = [constrViolation(maxIndex,:); constrViolation(r,:)];
    else
        constrViolationParents = [];
    end
    indicesParents = [maxIndex r];

function [indicesParents fitnessParents functionValueParents constrViolationParents equalityConstrValsParents inequalityConstrValsParents parents]=selectParents2(pop, fitnessVals, functionValue, constrViolation, equalityConstrVals, inequalityConstrVals, numParents)

    popSize = size(pop, 1);
    numOtherParents = numParents-1;
    
    [eliteFitness index]=max(fitnessVals);
    
    if ~isempty(constrViolation)
        t = find(constrViolation==0);
        if size(t,1)>2*numOtherParents
            permut = randperm(size(t,1));
            s = t(permut(1:2*numOtherParents));
        else
            [sortedFitness, I] = sort(fitnessVals, 'descend');
            s = I(1:2*numOtherParents);
        end
    else
        permut = randperm(popSize);
        s = permut(1:2*numOtherParents);
    end
    
    for i=1:2:2*numOtherParents
        if fitnessVals(s(i))>fitnessVals(s(i+1))
            r((i+1)/2) = s(i);
        else
            r((i+1)/2) = s(i+1);
        end
    end
        
    parents = [pop(index,:); pop(r,:)];
    fitnessParents = [fitnessVals(index); fitnessVals(r)];
    functionValueParents = [functionValue(index); functionValue(r)];

    if ~isempty(equalityConstrVals)
        equalityConstrValsParents = [equalityConstrVals(index,:); equalityConstrVals(r,:)];
    else
        equalityConstrValsParents = [];
    end
    if ~isempty(inequalityConstrVals)
        inequalityConstrValsParents = [inequalityConstrVals(index,:); inequalityConstrVals(r,:)];
    else
        inequalityConstrValsParents = [];
    end
    if ~isempty(constrViolation)
        constrViolationParents = [constrViolation(index,:); constrViolation(r,:)];
    else
        constrViolationParents = [];
    end
    indicesParents = [index r];
    
function offsprings = planar_pcx_crossover(parents, probCrossover, numOffsprings, param)

    %Planar PCX Crossover
    %Enter parents as a matrix of size [dim X numParents]
    %The code has been written such that the first parent is chosen as the
    %index parent. Any number of parents can be given as input. However,
    %the other two parents are chosen randomly. Only three parents can be
    %used to perform a crossover.

    sigma1 = .1;
    sigma2 = .1;
    dim = size(parents,2);
    numParents = size(parents,1);
    
    g = mean(parents);
    indexParent = parents(1,:);
    otherParents = parents(2:end,:);
    for i=1:numOffsprings
        otherParents = otherParents(randperm(numParents-1),:);
        if rand<probCrossover
            offsprings(i,:) = indexParent(1,:) + randn*sigma1*(indexParent(1,:) - g) + randn*sigma2*(otherParents(2,:) - otherParents(1,:))/2;
        else
            offsprings(i,:) = indexParent(1,:);
        end
    end
    
function offsprings = pcx_crossover(parents, probCrossover, numOffsprings, param)

    %PCX Crossover
    %Enter parents as a matrix of size [dim X numParents]
    %The code has been written such that the first parent is chosen as the
    %index parent. Any number of parents can be given as input.

    sigma1 = 0.1;
    sigma2 = 0.1;
    dim = size(parents,2);
    numParents = size(parents,1);
    
    g = mean(parents);

    d_g = g - parents(1,:);
    d_g_sq = abs(d_g).^2;
    mod_d_g = sum(d_g_sq,2).^(1/2);
    
    d = parents(2:numParents,:) - parents(ones(1,numParents-1),:);
    d_sq = abs(d).^2;
    mod_d = sum(d_sq,2).^(1/2);
    
    if mod_d_g == 0 | all(mod_d) == 0
        for i=1:numOffsprings
            offsprings(i,:) = parents(1,:);
        end
    else
        theta = dot(d,d_g(ones(1,numParents-1),:),2)./(mod_d.*mod_d_g);
        theta(theta>1) = 1;
        rho = mean(mod_d.*sqrt(1-theta.^2));
        for i=1:numOffsprings
            rand_vect = rho*randn(1,dim)*sigma1;
            mod_rand_vect = norm(rand_vect);

            perp_vect = rand_vect - d_g*dot(rand_vect,d_g)/(mod_d_g*mod_d_g);
            para_vect = randn*sigma2*d_g;
            offsprings(i,:) = parents(1,:) + perp_vect + para_vect;
        end
    end

function offsprings = ucx_crossover(parents, probCrossover, numOffsprings, param)

    %Uniform Centric Crossover
    %Enter parents as a matrix of size [dim X numParents]
    %The code has been written such that the first parent is chosen as the
    %index parent. Any number of parents can be given as input.
    alpha = param.alpha;
    kappa = param.kappa;
    constrViolation = param.constrViolation;
    popSize = param.popSize;
    
    if ~isempty(constrViolation)
        %sigma1 = sum(constrViolation==0)/popSize*alpha/2 + 0.2;
        sigma1 = sum(constrViolation==0)/popSize*alpha^.1 + 0.1;
    else
        %sigma1 = alpha/2 + 0.2;
        sigma1 = alpha^.1 + 0.1;
    end
    %sigma1 = sum(constrViolation==0)/popSize*(kappa^(-(1-alpha)) - 1/kappa) + .2; %Another strategy
    %sigma2 = sigma1/2;
    sigma2 = sigma1;
    dim = size(parents,2);
    numParents = size(parents,1);
    
    g = mean(parents(1:end,:));

    d_g = g - parents(1,:);
    d_g_sq = abs(d_g).^2;
    mod_d_g = sum(d_g_sq,2).^(1/2);
    
    d = parents(2:numParents,:) - parents(ones(1,numParents-1),:);
    d_sq = abs(d).^2;
    mod_d = sum(d_sq,2).^(1/2);
    
    if mod_d_g == 0 | all(mod_d) == 0
        for i=1:numOffsprings
            offsprings(i,:) = mutation(parents(1,:),1);
        end
    else
        theta = dot(d,d_g(ones(1,numParents-1),:),2)./(mod_d.*mod_d_g);
        theta(theta>1) = 1;
        rho = mean(mod_d.*sqrt(1-theta.^2));
        for i=1:numOffsprings
            rand_vect = rho*randn(1,dim)*sigma1;
            mod_rand_vect = norm(rand_vect);
            perp_vect = rand_vect - d_g*dot(rand_vect,d_g)/(mod_d_g*mod_d_g);
            para_vect = randn*sigma2*d_g;
            offsprings(i,:) = parents(1,:) + (2*rand)*(g-parents(1,:)) + perp_vect + para_vect;
        end
    end

function offsprings = adaptive_planar_pcx_crossover(parents, probCrossover, numOffsprings, param)

    %Planar PCX Crossover
    %Enter parents as a matrix of size [dim X numParents]
    %The code has been written such that the first parent is chosen as the
    %index parent. Any number of parents can be given as input. However,
    %the other two parents are chosen randomly. Only three parents can be
    %used to perform a crossover.
    epsilon = 1e-10;
    dim = size(parents,2);
    numParents = size(parents,1);
    
    g = mean(parents);
    indexParent = parents(1,:);
    otherParents = parents(2:end,:);
    
    sigma1 = .1;
    if mean(abs(indexParent(1,:) - g)) <= epsilon
        sigma2 = 0;
    else
        sigma2 = 1/mean(abs(indexParent(1,:) - g));
    end
    
    for i=1:numOffsprings
        otherParents = otherParents(randperm(numParents-1),:);
        if rand<probCrossover
            offsprings(i,:) = indexParent(1,:) + randn*sigma1*(indexParent(1,:) - g) + randn*sigma2*(otherParents(2,:) - otherParents(1,:))/2;
        else
            offsprings(i,:) = indexParent(1,:);
        end
    end
    
function offsprings = mutation(offsprings, probMutation)

    numOffsprings=size(offsprings,1);
    dim=size(offsprings,2);
    mum=20;
    for i=1:numOffsprings
        r = rand(1,dim);
        index = r<0.5;
        delta(index) = (2*r(index)).^(1/(mum+1)) - 1;
        index = r>=0.5;
        delta(index) = 1 - (2*(1 - r(index))).^(1/(mum+1));

        % Generate the corresponding child element.
        r = rand(1,dim);
        if ~all(r>=probMutation)
            offsprings(i,r<probMutation) = offsprings(i,r<probMutation) + delta(r<probMutation);
        end
    end
    
function offsprings = mutationVarianceBased(offsprings, probMutation, variance)

    %Variance based mutation
    numOffsprings=size(offsprings,1);
    dim=size(offsprings,2);
    mum=20;
    if isempty(variance)
        variance = zeros(1,dim);
    end
    for i=1:numOffsprings
        r = rand(1,dim);
        index = r<0.5;
        delta(index) = (2*r(index)).^(1/(mum+1)) - 1;
        index = r>=0.5;
        delta(index) = 1 - (2*(1 - r(index))).^(1/(mum+1));

        % Generate the corresponding child element.
        r = rand(1,dim);
        if ~all(r>=probMutation)
            offsprings(i,r<probMutation) = offsprings(i,r<probMutation) + (1+sqrt(variance(r<probMutation))).*delta(r<probMutation);
        end
    end
    
function [pop fitnessVals functionValue constrViolation equalityConstrVals inequalityConstrVals] = update(pop, fitnessVals, functionValue, constrViolation, equalityConstrVals, inequalityConstrVals, indicesParents, parents, fitnessParents, functionValueParents, constrViolationParents, equalityConstrValsParents, inequalityConstrValsParents, offsprings, fitnessOffsprings, functionValueOffsprings, constrViolationOffsprings,  equalityConstrValsOffsprings, inequalityConstrValsOffsprings)

    chooseMembers = 2;
    numOffsprings = size(offsprings,1);
    popSize = size(pop,1);
    permut = randperm(popSize);
    r = permut(1:chooseMembers);
    
    pool = [pop(r,:); offsprings];
    fitnessPool = [fitnessVals(r); fitnessOffsprings];
    functionValuePool = [functionValue(r); functionValueOffsprings];
    if ~isempty(constrViolation) & ~isempty(constrViolationOffsprings)
        constrViolationPool = [constrViolation(r); constrViolationOffsprings];
        if ~isempty(equalityConstrVals) && ~isempty(equalityConstrValsOffsprings)
            equalityConstrValsPool = [equalityConstrVals(r,:); equalityConstrValsOffsprings];
        else
            equalityConstrValsPool = [];
        end
        if ~isempty(inequalityConstrVals) && ~isempty(inequalityConstrValsOffsprings)
            inequalityConstrValsPool = [inequalityConstrVals(r,:); inequalityConstrValsOffsprings];
        else
            inequalityConstrValsPool = [];
        end
    else
        constrViolationPool = [];
    end
    
    [tempFitnessPool,I] = sort(fitnessPool,'descend');
    
    pop(r,:) = pool(I(1:chooseMembers),:);
    fitnessVals(r) = fitnessPool(I(1:chooseMembers));
    functionValue(r) = functionValuePool(I(1:chooseMembers));
    if ~isempty(constrViolationPool)
        constrViolation(r) = constrViolationPool(I(1:chooseMembers));
        if ~isempty(equalityConstrVals)
            equalityConstrVals(r,:) = equalityConstrValsPool(I(1:chooseMembers),:);
        end
        if ~isempty(inequalityConstrVals)
            inequalityConstrVals(r,:) = inequalityConstrValsPool(I(1:chooseMembers),:);
        end
    end

function offsprings=checkLimits(offsprings, dimMin, dimMax)

    numOffsprings = size(offsprings,1);
    dimMinMatrix = dimMin(ones(1,numOffsprings),:);
    offsprings(offsprings<dimMinMatrix)=dimMinMatrix(offsprings<dimMinMatrix);
    dimMaxMatrix = dimMax(ones(1,numOffsprings),:);
    offsprings(offsprings>dimMaxMatrix)=dimMaxMatrix(offsprings>dimMaxMatrix);
    
function offsprings=checkLimitsReflection(offsprings, dimMin, dimMax)

    %This function reflects an infeasible point into the variable bounds. If the
    %point lies far away, it assigns it a random position in the bounds.
    numOffsprings = size(offsprings,1);
    dimMinMatrix = dimMin(ones(1,numOffsprings),:);
    dimMaxMatrix = dimMax(ones(1,numOffsprings),:);
    i = 0;
    while sum(sum(offsprings<dimMinMatrix)) | sum(sum(offsprings>dimMaxMatrix))
        I = offsprings<dimMinMatrix-(dimMaxMatrix-dimMinMatrix);
        J = offsprings>dimMaxMatrix+(dimMaxMatrix-dimMinMatrix);
        randI = rand(size(I));
        randJ = rand(size(J));
        offsprings(I) = dimMinMatrix(I) + randI(I).*(dimMaxMatrix(I)-dimMinMatrix(I));
        offsprings(J) = dimMinMatrix(J) + randJ(J).*(dimMaxMatrix(J)-dimMinMatrix(J));
        offsprings(offsprings<dimMinMatrix)=offsprings(offsprings<dimMinMatrix) + 2*(dimMinMatrix(offsprings<dimMinMatrix)-offsprings(offsprings<dimMinMatrix));
        offsprings(offsprings>dimMaxMatrix)=offsprings(offsprings>dimMaxMatrix) + 2*(dimMaxMatrix(offsprings>dimMaxMatrix)-offsprings(offsprings>dimMaxMatrix));
    end
    
function [fitnessVals constrViolation] = assignFitness(functionValue, equalityConstrVals, inequalityConstrVals, violationWindowEquality, violationWindowInequality)
 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Specify the number of constraints
    numEqualityConstr = size(equalityConstrVals,2);
    numInequalityConstr = size(inequalityConstrVals,2);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if numEqualityConstr == 0 & numInequalityConstr == 0
        constrViolation = [];
    end
    if numEqualityConstr == 0 & numInequalityConstr ~= 0
        inequalityConstrVals = inequalityConstrVals - violationWindowInequality;
        inequalityConstrVals(inequalityConstrVals<0)=0;
        constrViolation = sum(inequalityConstrVals,2);
    end
    if numEqualityConstr ~= 0 & numInequalityConstr == 0
        equalityConstrVals(abs(equalityConstrVals)<=violationWindowEquality) = 0;
        equalityConstrVals(abs(equalityConstrVals)>violationWindowEquality) = equalityConstrVals(abs(equalityConstrVals)>violationWindowEquality)-violationWindowEquality;
        constrViolation = sum(abs(equalityConstrVals),2);
    end
    if numEqualityConstr ~= 0 & numInequalityConstr ~= 0
        inequalityConstrVals = inequalityConstrVals - violationWindowInequality;
        inequalityConstrVals(inequalityConstrVals<0)=0;
        
        equalityConstrVals(abs(equalityConstrVals)<=violationWindowEquality) = 0;
        equalityConstrVals(abs(equalityConstrVals)>violationWindowEquality) = equalityConstrVals(abs(equalityConstrVals)>violationWindowEquality)-violationWindowEquality;
        constrViolation = sum(inequalityConstrVals,2)+ sum(abs(equalityConstrVals),2);
    end

    if (isempty(constrViolation))
        fitnessVals = functionValue;
        return;
    end
    
    if (constrViolation>0)
        fitnessVals = -constrViolation;
        return;
    end
    
    %Otherwise
    fitnessVals(constrViolation>0,1) = min(functionValue(constrViolation==0)) - constrViolation(constrViolation>0);
    fitnessVals(constrViolation==0,1) = functionValue(constrViolation==0);
    
function approxFunctionValue = approximatedFunction(pop, parameters)

    approxFunctionValue = parameters.constant + pop*parameters.linear + pop*parameters.sqmatrix*pop';

function [c ceq] = approximatedConstraints(pop, parametersEqualityConstr, parametersInequalityConstr)

    if length(parametersEqualityConstr)~=0
        for i=1:length(parametersEqualityConstr)
            ceq(i) = parametersEqualityConstr{i}.constant + pop*parametersEqualityConstr{i}.linear + pop*parametersEqualityConstr{i}.sqmatrix*pop';
        end
    else
        ceq = [];
    end
    
    if length(parametersInequalityConstr)~=0
        for i=1:length(parametersInequalityConstr)
            c(i) = parametersInequalityConstr{i}.constant + pop*parametersInequalityConstr{i}.linear + pop*parametersInequalityConstr{i}.sqmatrix*pop';
        end
    else
        c = []; %Negative value suggests that there is no active inequality
    end
    
function bool = isLinear(approximateFunctionParameters,epsilonZero)

    bool = true;
    if sum(sum(abs(approximateFunctionParameters.function.sqmatrix)>epsilonZero))>0
        bool = false;
    end
    for i=1:length(approximateFunctionParameters.equalityConstr)
        if sum(sum(abs(approximateFunctionParameters.equalityConstr{i}.sqmatrix)>epsilonZero))>0
            bool = false;
        end
    end
    for i=1:length(approximateFunctionParameters.inequalityConstr)
        if sum(sum(abs(approximateFunctionParameters.inequalityConstr{i}.sqmatrix)>epsilonZero))>0
            bool = false;
        end
    end
    
function bool = isFeasible(member,constrViolation,dimMin,dimMax,violationWindowConstraint)

    bool = 0;
    if constrViolation>violationWindowConstraint
        return;
    end
    if member>=dimMin
        if member<=dimMax
            bool = 1;
            return;
        end
    end