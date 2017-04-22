function [eliteFunctionValue, llEliteFunctionValue, eliteIndiv, llEliteIndiv, ulFunctionEvaluations, llFunctionEvaluations, llCalls, gen, stoppingCondition,outputLowerRunTags, outputLowerRunGens]=ulSearch(testProblemName, ulPopSize, ulMaxGens, ulDim, ulDimMin, ulDimMax, llPopSize, llMaxGens, llDim, llDimMin, llDimMax, ulEpsilonStopping, llEpsilonStopping, conditionCheck, seed)
 
%Call as [eliteFitness, eliteIndiv]=g3pcx(testProblemName, popSize,
%maxGens, dim, dimMin, dimMax)
%Algorithm for maximization problem
%feasibility g(x)<0 and h(x)=0
 
warning off all;
global ulFunctionEvaluations 
global llFunctionEvaluations;

ulFunctionEvaluations = 0;
llFunctionEvaluations = 0;

if nargin==15
    rng(seed);
end

if nargin>=14
    if strcmp(conditionCheck,'continue') == 1
        previousRun = load([testProblemName '_temp.mat']);
    end
end

if nargin==11
    ulEpsilonStopping = 1e-4;
    llEpsilonStopping = 1e-5;
end

if nargin<11
    display('Insufficient input arguments. Exiting.')
    return;
end

stoppingCondition = [];
archive = struct('tag1',[],'tag0',[]); % Tag 1 member Archive used for meta-model construction

probCrossover=0.9;         % The probability of crossover 
probMutation=1/ulDim;      % The mutation probability
%probMutation=0.1;          % The mutation probability

parentSelectionType = 2;   % 1 => selectParents
                           % 2 => selectParents2
                           % 3 => tournamentSelectParents
                           % 4 => tournamentSelectSortedParents

crossoverType=3;           % 0 => no crossover
                           % 1 => planar pcx crossover
                           % 2 => normal pcx crossover
                           % 3 => adaptive planar pcx crossover
 
mutationType=1;            % 0 => no mutation
                           % 1 => polynomial mutation
                           % 2 => variance based polynomial mutation
 
visualizationFlag=2;       % 0 => don't visualize
                           % 1 => visualize for all
                           % 2 => visualize only for multiobjective problems
 
verboseFlag=1;             % 0 => run quietly
                           % 1 => display details of each generation
 
maxError = 1e-4;
minTagOne = ceil(0.5*ulPopSize);

numParents=3;
numOffsprings=2;
llLocalSearchPopSize = max([4 llPopSize/10]);
llLocalSearchMaxGens = llMaxGens;
 
probSearch=1.0;            % Lower level search probability
llCalls = 0;
violationWindowEqualityInitial = 1e-2;
violationWindowInequalityInitial = 1e-4;
 
eliteIndiv=[];
eliteFitness=-realmax;

minPopSize = (ulDim+1)*(ulDim+2)/2+5*(ulDim);

if ulPopSize < minPopSize
    display(['Population size at the upper level was too small, automatically increased by the code to ' minPopSize])
    ulPopSize = minPopSize;
end
 
if nargin>=14
    if strcmp(conditionCheck,'continue') == 1
        ulPop=initialization(ulDim, ulDimMin, ulDimMax, llDimMin, llDimMax, ulPopSize-1, testProblemName);
        ulPop(end+1,:) = previousRun.ulPop(previousRun.maxIndex,:);
    end
end
        
if (nargin<14) || (strcmp(conditionCheck,'continue') == 0)
    ulPop = initialization(ulDim, ulDimMin, ulDimMax, llDimMin, llDimMax, ulPopSize, testProblemName);
    %ulPop =ulInitialize(ulPopSize,ulDim,ulDimMin,ulDimMax);
end

% initialize popultion cell array 
population = cell(ulPopSize,1);
ulInitializationRecheck = 0; % do not re-evaluate the initial population

for i=1:ulPopSize
    display(['Doing lower level search for population member ' num2str(i)]);
    [llPop(i,:), llFunctionValue, llConstr, tag.ulPop(i), ~, lowerGen]=llSearch(testProblemName, llPopSize, llMaxGens, llDim, llDimMin, llDimMax, ulPop(i,:), llEpsilonStopping, [], [], []);
    %Note upper level tag 1 member ensures only lower level optimality, the
    %member may still be infeasible because of lower level constraints
    
    llCalls = llCalls+1;
    outputLowerRunTags(llCalls) = tag.ulPop(i);
    outputLowerRunGens(llCalls) = lowerGen;
    if lowerGen == 0
        optType = 'Quadratic Programming';
    else
        optType = 'Quadratic Programming + GA';
        ulInitializationRecheck = 1;
    end
    
    % upper level function and constraint value calculation
    [functionValue, equalityConstrVals, inequalityConstrVals]=ulTestProblem(ulPop(i,:), llPop(i,:), testProblemName);
    
    population{i} = struct('upper',ulPop(i,:),'lower',llPop(i,:),'functionValue',functionValue,'llFunctionValue',llFunctionValue, ...
                                    'fitness',[],'constrViolation',[],'equalityConstrVals',equalityConstrVals,'inequalityConstrVals',inequalityConstrVals,...
                                    'llEqualityConstrVals',llConstr.equalityConstrVals,'llInequalityConstrVals',llConstr.inequalityConstrVals,'optType',optType);  
end

% re-evaluate the initial population (coevolution)
[~,centroids] = kmeans(llPop,llLocalSearchPopSize,'EmptyAction','singleton');
if ulInitializationRecheck == 1   
    for i=1:ulPopSize
        display(['Doing lower level search for population member ' num2str(i)]);
        [llPop(i,:), llFunctionValue, llConstr, tag.ulPop(i), ~, lowerGen]=llSearch(testProblemName, llLocalSearchPopSize, llLocalSearchMaxGens, llDim, llDimMin, llDimMax, ulPop(i,:), llEpsilonStopping, [], [], centroids);
        %Note upper level tag 1 member ensures only lower level optimality, the
        %member may still be infeasible because of lower level constraints

        llCalls = llCalls+1;
        outputLowerRunTags(llCalls) = tag.ulPop(i);
        outputLowerRunGens(llCalls) = lowerGen;
        if lowerGen == 0
            optType = 'Quadratic Programming';
        else
            optType = 'Quadratic Programming + GA';
        end

        % upper level function and constraint value calculation
        [functionValue, equalityConstrVals, inequalityConstrVals]=ulTestProblem(ulPop(i,:), llPop(i,:), testProblemName);

        population{i} = struct('upper',ulPop(i,:),'lower',llPop(i,:),'functionValue',functionValue,'llFunctionValue',llFunctionValue, ...
                                        'fitness',[],'constrViolation',[],'equalityConstrVals',equalityConstrVals,'inequalityConstrVals',inequalityConstrVals,...
                                        'llEqualityConstrVals',llConstr.equalityConstrVals,'llInequalityConstrVals',llConstr.inequalityConstrVals,'optType',optType);  
    end
end

%Minimum number of lower level points required for quadratic approximation
%of the optimal lower level values
numberLowerLevelPointsTrain = (ulDim+1)*(ulDim+2)/2+2*(ulDim);
numberLowerLevelPointsEval = ulDim;
minSizeMappings = numberLowerLevelPointsTrain + numberLowerLevelPointsEval;
minSizePhiLS = (ulDim+llDim+1)*(ulDim+llDim+2)/2 + 2*(ulDim) + ulDim;

archiveSize = 10*minSizeMappings;

llOptimaFitness(1:ulPopSize,:) = 0;
llOptimaFitnessOffsprings(1:numOffsprings,:) = 0;
gen = 0;

% evaluate the fitness of the population. The vector of fitness values 
% returned  must be of dimensions ulPopSize x 1.
[~, ~, population] = assignFitness(population, violationWindowEqualityInitial, violationWindowInequalityInitial, tag.ulPop);

archive.tag1 = archiveUpdate(archive.tag1,population(tag.ulPop == 1),archiveSize);

alphaStoppingInitial = sum(var([ulPop(tag.ulPop==1,:) llPop(tag.ulPop==1,:)]))/(ulDim+llDim);
llMemberVariance = var(llPop); 
% ulMemberVariance = var(ulPop);
 
for gen=1:ulMaxGens
    % reduce equality and inequality 
    violationWindowEquality = violationWindowEqualityInitial*(1+exp(-gen*numOffsprings/ulPopSize));
    violationWindowInequality = violationWindowInequalityInitial*(1+exp(-gen*numOffsprings/ulPopSize));
    % select parents from current population to performance recombination
    [indicesParents, ~, ~, ~, ~, ~, parents] = parentSelection(population, parentSelectionType, numParents, tag);
	tag.parents = tag.ulPop(indicesParents);
    
    % perform recombination (crossover + mutation)
    offsprings = geneticOperation (parents,probCrossover,probMutation,numOffsprings,crossoverType,mutationType);
    offsprings = checkLimits(offsprings, ulDimMin, ulDimMax);
    
    parentsLowerLevelVariables = llPop(indicesParents,:);
    poDistance = computeDistance(offsprings, parents);
    children = cell(numOffsprings,1);
    
    for i=1:numOffsprings
        % lower level optimization
        if rand<probSearch || sum(tag.ulPop==1)<minTagOne || length(archive.tag1)<minSizeMappings
            [~, closestParent] = sort(poDistance(i,:),'ascend');
            taggedClosestParent = closestParent(1);
            for j=1:length(closestParent)
                if tag.parents(closestParent(j))==1
                    taggedClosestParent = closestParent(j);
                    break;
                end
            end
            offspringsLowerLevelVariables = parentsLowerLevelVariables(taggedClosestParent,:);
            [offspringsLowerLevelVariables, llFunctionValue, llConstr, tag.offsprings(i), ~, lowerGen]=llSearch(testProblemName, llLocalSearchPopSize, llLocalSearchMaxGens, llDim, llDimMin, llDimMax, offsprings(i,:), llEpsilonStopping, offspringsLowerLevelVariables, llMemberVariance, llPop);
            llEqualityConstrVals = llConstr.equalityConstrVals; 
            llInequalityConstrVals = llConstr.inequalityConstrVals;
            if tag.offsprings(i)==0
            	llOptimaFitnessOffsprings(i) = -realmax;
            else
            	llOptimaFitnessOffsprings(i) = 0;
            end
            llCalls = llCalls+1;
            outputLowerRunTags(llCalls) = tag.offsprings(i);
            outputLowerRunGens(llCalls) = lowerGen;
            if lowerGen == 0
            	optTypeOffsprings = 'Quadratic Programming';
            else
            	optTypeOffsprings = 'Quadratic Programming + GA';
            end
        else
            % offsprings lower level optimal variables using mapping
            % approximation
            [psiMapping,phiMapping,lies] = getMappings(offsprings(i,:),archive.tag1);
            [offspringsLowerLevelVariables,optTypeOffsprings,sumMSE,validMSE] = getLowerLevelVariableFromMapping(offsprings(i,:),psiMapping,phiMapping,ulDim,llDim,archive);
            llOptimaFitnessOffsprings(i) = -sumMSE;
            %If lower level variables lie outside the variable bounds it
            %should be tagged 0
            lowerLevelLimitCheck = 0;
            if (offspringsLowerLevelVariables-checkLimits(offspringsLowerLevelVariables, llDimMin, llDimMax))==0
                lowerLevelLimitCheck = 1;
            end
            offspringsLowerLevelVariables=checkLimits(offspringsLowerLevelVariables, llDimMin, llDimMax);
            if (lies==1 && lowerLevelLimitCheck==1 && validMSE<maxError)
                tag.offsprings(i) = 1;
            else
            	tag.offsprings(i) = 0;
            end
            [llFunctionValue,llEqualityConstrVals,llInequalityConstrVals] = llTestProblem(offspringsLowerLevelVariables,testProblemName,offsprings(i,:));  
        end
          
        [functionValueOffsprings, equalityConstrValsOffsprings, inequalityConstrValsOffsprings]=ulTestProblem(offsprings(i,:), offspringsLowerLevelVariables, testProblemName);
        
        children{i} = struct('upper',offsprings(i,:),'lower',offspringsLowerLevelVariables,'functionValue',functionValueOffsprings,'llFunctionValue',llFunctionValue, ...
                                    'fitness',[],'constrViolation',[],'equalityConstrVals',equalityConstrValsOffsprings,'inequalityConstrVals',inequalityConstrValsOffsprings,...
                                    'llEqualityConstrVals',llEqualityConstrVals,'llInequalityConstrVals',llInequalityConstrVals,'optType',optTypeOffsprings);
    end
    
    [~, ~, tempPopulation] = assignFitness([population;children],violationWindowEquality, violationWindowInequality, [tag.ulPop tag.offsprings]);
    population = tempPopulation(1:ulPopSize);
    children = tempPopulation(1+ulPopSize:end);
    
    % update archive with children depending on tag 
	if sum(tag.offsprings == 1) > 0 %#ok<ALIGN>
        archive.tag1 = archiveUpdate(archive.tag1,children(tag.offsprings == 1),archiveSize);
    end
    if sum(tag.offsprings == 0) > 0
    	archive.tag0 = archiveUpdate(archive.tag0,children(tag.offsprings == 0),archiveSize);
    end
    
    [population,fitnessVals,llOptimaFitness, tag] = update(population,children,llOptimaFitness, llOptimaFitnessOffsprings,tag,tag.offsprings);
    llPop = cell2mat(cellfun(@(x) x.lower, population,'UniformOutput',false)); llMemberVariance = var(llPop);
    ulPop = cell2mat(cellfun(@(x) x.upper, population,'UniformOutput',false)); ulMemberVariance = var(ulPop);
    
    % identifying the elite member 
    if sum(tag.ulPop==1)==0
        [eliteFitness, maxIndex]=max(fitnessVals);
        avgFitness = mean(fitnessVals);
    else
        I = find(tag.ulPop==1);
        [eliteFitness, index]=max(fitnessVals(tag.ulPop==1));
        maxIndex = I(index);
        avgFitness = mean(fitnessVals(tag.ulPop == 1));
    end
    eliteMember = population{maxIndex}; % elite member cell of struct
    eliteIndiv=eliteMember.upper;
    llEliteIndiv=eliteMember.lower;
    eliteFunctionValue=eliteMember.functionValue;
    llEliteFunctionValue = eliteMember.llFunctionValue;
    if ~isempty(eliteMember.constrViolation)
        eliteConstrViolation = eliteMember.constrViolation;
    else
        eliteConstrViolation=0;
    end
    
    % display the generation number, the average Fitness of the population,
    % and the maximum function value in the population
    if verboseFlag
        if (eliteConstrViolation>0)
            display(['gen=' num2str(gen,'%.3d') '   No feasible solution found.']);
        elseif sum(tag.ulPop==1)==0
            display(['gen=' num2str(gen,'%.3d') '   No feasible solution found.']);
        else
            if size(functionValue,2)==1
                display(['gen=' num2str(gen,'%.3d') '   avgFitness=' ...
                num2str(avgFitness,'%3.3f') '   eliteFitnessValue=' ...
                num2str(eliteFitness,'%3.3f')]);
                display(['eliteFunctionValue=' num2str(eliteFunctionValue)]);
            else
                display('Bilevel feasible solutions found. Improving the non-dominated set.')
            end
        end
    end
    % Perform visualization
    if visualizationFlag>0 && mod(gen,10)==0
       if size(functionValue,2)==1 && visualizationFlag==1
            figure(1)
            set (gcf, 'color', 'w');
            hold off
            plot(1:ulPopSize,fitnessVals, '.');
            title(['Generation = ' num2str(gen) ', Average Fitness = ' sprintf('%0.3f', avgFitness)]);
            ylabel('Fitness');
            xlabel('Population members');
       end
        if size(functionValue,2)>=2 && visualizationFlag>=2
            figure(1)
            set (gcf, 'color', 'w');
            hold off
            functionValueForPlot = cell2mat(cellfun(@(x) x.functionValue, population,'UniformOutput',false));
            plot(functionValueForPlot(:,1),functionValueForPlot(:,2),'*');
            title(['Population at Generation = ' num2str(gen)]);
            xlabel('Objective 1');
            ylabel('Objective 2');
        end
        drawnow;
    end
    
    eliteFunctionValueAtGen(gen,:) = eliteFunctionValue;
    stoppingParameters.ulEpsilonStopping = ulEpsilonStopping;
    stoppingParameters.llEpsilonStopping = llEpsilonStopping;
    stoppingParameters.alphaStoppingInitial = alphaStoppingInitial;
    stoppingParameters.eliteFunctionValueAtGen = eliteFunctionValueAtGen;
    stoppingParameters.llEliteFunctionValue = llEliteFunctionValue;
    stoppingParameters.eliteConstrViolation = eliteConstrViolation; 
    stoppingParameters.eliteIndiv = eliteIndiv;
    stoppingParameters.llEliteIndiv = llEliteIndiv;
    stoppingParameters.testProblemName = testProblemName;
    [StoppingCriteria, stoppingCondition] = terminationCheck(gen, tag, ulPop, llPop, ulDim, llDim, stoppingParameters);
    
    % improvements - local search + recheck
    Flag = improvementsFlag(ulPopSize,gen,ulMaxGens,StoppingCriteria,archive,minSizeMappings,minSizePhiLS);
    %Do not do local search if the problem is multi-objective at upper
    %level
    if size(eliteFunctionValue,2)>1
        Flag.localSearch = 0;
    end
    if Flag.localSearch > 0
        if rand < 0.5
            initialIndivLS = eliteMember;
        else
            initialIndivLS = children{randi(numOffsprings)};
        end
        if exist('localSearchTemp','var')
            if strcmp(localSearchTemp.method,'Approx') || strcmp(localSearchTemp.method,'ApproxPhi')
                if localSearchTemp.tagAcceptLS == 0
                    Flag.localSearch = 1;
                end
            end
        end
     
        [eliteIndivLS, llEliteIndivLS,localSearchTemp] = doLocalSearch2(archive, initialIndivLS,ulDimMin, ulDimMax, llDimMin, llDimMax,Flag,testProblemName);
        [llEliteIndivLS, llEliteFunctionValueLS, llEliteConstrLS, tag.lowerIndivLS, ~, lowerGenLS]=llSearch(testProblemName, llLocalSearchPopSize, llLocalSearchMaxGens, llDim, llDimMin, llDimMax, eliteIndivLS, llEpsilonStopping, llEliteIndivLS,llMemberVariance,llPop);
        llCalls = llCalls+1;
        outputLowerRunTags(llCalls) = tag.lowerIndivLS;
        outputLowerRunGens(llCalls) = lowerGenLS;
        [eliteFunctionValueLS, eliteEqualityConstrValsLS, eliteInequalityConstrValsLS]=ulTestProblem(eliteIndivLS, llEliteIndivLS, testProblemName);
        if lowerGenLS == 0
            optTypeLS = 'Quadratic Programming';
        else
            optTypeLS = 'Quadratic Programming + GA';
        end
        eliteLS = struct('upper',eliteIndivLS,'lower',llEliteIndivLS,'functionValue',eliteFunctionValueLS,'llFunctionValue',llEliteFunctionValueLS, ...
                                    'fitness',[],'constrViolation',[],'equalityConstrVals',eliteEqualityConstrValsLS,'inequalityConstrVals',eliteInequalityConstrValsLS,...
                                    'llEqualityConstrVals',llEliteConstrLS.equalityConstrVals,'llInequalityConstrVals',llEliteConstrLS.inequalityConstrVals,'optType',optTypeLS);
        if tag.lowerIndivLS == 1
            Flag.update = 1;
            llOptimaFitnessLS = 0;
            localSearchTemp.tagAcceptLS = ifAcceptLS(initialIndivLS, eliteLS, violationWindowEquality, violationWindowInequality);
        else
            localSearchTemp.tagAcceptLS = 0;
        end
    end
    
    if Flag.recheck == 1
        [llPopCheck, llFunctionValueCheck, llConstrCheck, tag.ulPopCheck, ~, lowerGenCheck]=llSearch(testProblemName,  llPopSize, llMaxGens, llDim, llDimMin, llDimMax, eliteIndiv, llEpsilonStopping, [], [], []);
        llCalls = llCalls+1;
        outputLowerRunTags(llCalls) = tag.ulPopCheck;
        outputLowerRunGens(llCalls) = lowerGen;
        if lowerGenCheck == 0
            optTypeCheck = 'Quadratic Programming';
        else
            optTypeCheck = 'Quadratic Programming + GA';
        end
        [functionValueCheck, equalityConstrValsCheck, inequalityConstrValsCheck]=ulTestProblem(eliteIndiv, llPopCheck, testProblemName);
        
        eliteReplacement = struct('upper',eliteIndiv,'lower',llPopCheck,'functionValue',functionValueCheck,'llFunctionValue',llFunctionValueCheck, ...
                                    'fitness',[],'constrViolation',[],'equalityConstrVals',equalityConstrValsCheck,'inequalityConstrVals',inequalityConstrValsCheck,...
                                    'llEqualityConstrVals',llConstrCheck.equalityConstrVals,'llInequalityConstrVals',llConstrCheck.inequalityConstrVals,'optType',optTypeCheck);  
        if llEliteFunctionValue<llFunctionValueCheck && tag.ulPopCheck==1
            Flag.replace = 1;
        end
    end
    
    % improvements - corrections if needed 
    if Flag.update == 1
        [~, ~, tempPopulation] = assignFitness([population; eliteLS], violationWindowEquality, violationWindowInequality*(1+exp(-gen*numOffsprings/ulPopSize)));
        population = tempPopulation(1:ulPopSize);
        eliteLS = tempPopulation(1+ulPopSize:end);
        % include the LS point into archive
        archive.tag1 = archiveUpdate(archive.tag1,eliteLS,archiveSize,1);
        % update the current population
        [population,fitnessVals,llOptimaFitness, tag] = update(population,eliteLS,llOptimaFitness, llOptimaFitnessLS,tag,tag.lowerIndivLS);
        llPop = cell2mat(cellfun(@(x) x.lower, population,'UniformOutput',false)); llMemberVariance = var(llPop);
        ulPop = cell2mat(cellfun(@(x) x.upper, population,'UniformOutput',false)); ulMemberVariance = var(ulPop);
    end
    
    if Flag.replace == 1
        % correct the population member 
        replaceIndex = maxIndex;
        tag.ulPop(replaceIndex) =  tag.ulPopCheck;
        population{replaceIndex} = eliteReplacement;
        eliteFunctionValueAtGen(gen,:) = eliteReplacement.functionValue;
        
%         %Makes all upper level members worse than current member 
%         tag.ulPop(functionValueCheck < cellfun(@(x) x.functionValue, population)) = 0;
%         archive(functionValueCheck < cellfun(@(x) x.functionValue, archive)) = [];
        [~, ~, population] = assignFitness(population, violationWindowEquality, violationWindowInequality, tag.ulPop);    
        % assuming the elite member in the current population is the same
        % as the elite member in the archive
        archive.tag1 = archiveUpdate(archive.tag1,population{replaceIndex},archiveSize);
        llPop = cell2mat(cellfun(@(x) x.lower, population,'UniformOutput',false)); llMemberVariance = var(llPop);
        ulPop = cell2mat(cellfun(@(x) x.upper, population,'UniformOutput',false)); ulMemberVariance = var(ulPop);
    end
    
    if StoppingCriteria==1
        disp([stoppingCondition,' termination criterion met']);
        % identifying the elite member 
        if sum(tag.ulPop==1)==0
            [~, maxIndex]=max(fitnessVals);
        else
            I = find(tag.ulPop==1);
            [~, index]=max(fitnessVals(tag.ulPop==1));
            maxIndex = I(index);
        end
        
        if size(functionValue,2)==1
            eliteMember = population{maxIndex}; % elite member cell of struct
            eliteIndiv=eliteMember.upper;
            llEliteIndiv=eliteMember.lower;
            eliteFunctionValue=eliteMember.functionValue;
            llEliteFunctionValue = eliteMember.llFunctionValue;
        else
            for i=1:length(population)
                eliteIndiv(i,:)=population{i}.upper;
                llEliteIndiv(i,:)=population{i}.lower;
                eliteFunctionValue(i,:)=population{i}.functionValue;
                llEliteFunctionValue(i,:)=population{i}.llFunctionValue;
            end
        end
        
        ulFunctionEvaluations
        llFunctionEvaluations
        if size(functionValue,2)==1
            display(['No of lower level calls = ',num2str(llCalls)]);
            display(['Upper level function value = ',num2str(eliteFunctionValue)]);
            display(['Lower level function value = ',num2str(llEliteFunctionValue)]);
            display(['Upper level elite vector = ',num2str(eliteIndiv)]);
            display(['Lower level elite vector = ',num2str(llEliteIndiv)]);
        else
            display('No of lower level calls = '); llCalls
            display('Upper level function values = '); eliteFunctionValue
            display('Lower level function values = '); llEliteFunctionValue
            display('Upper level elite vectors = '); eliteIndiv
            display('Lower level elite vectors = '); llEliteIndiv
        end
        return;
    end
    save([testProblemName '_temp'], 'ulPopSize', 'ulMaxGens', 'ulDim', 'ulDimMin', 'ulDimMax', 'llPopSize', 'llMaxGens', 'llDim', 'llDimMin', 'llDimMax', 'ulPop', 'llPop', 'maxIndex');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%End of main function%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ulPop=ulInitialize(ulPopSize,ulDim,ulDimMin,ulDimMax)
 
    for i=1:ulPopSize
        ulPop(i,:) = ulDimMin + rand(1, ulDim).*(ulDimMax-ulDimMin);
    end
    
function d = computeDistance(matrix1, matrix2)
 
    %Computes pairwise distance between rows of matrix1 and matrix2
    sz1 = size(matrix1, 1);
    sz2 = size(matrix2, 1);
    
    for i = 1:sz1
        for j = 1:sz2
            d(i,j) = sqrt(sum((matrix1(i,:)-matrix2(j,:)).^2));
        end
    end
    
function [indicesParents, fitnessParents, functionValueParents, constrViolationParents, equalityConstrValsParents, inequalityConstrValsParents, parents] = parentSelection(population, parentSelectionType, numParents, tag)
    
    ulPop = cell2mat(cellfun(@(x) x.upper, population,'UniformOutput',false));
    fitnessVals = cellfun(@(x) x.fitness, population);
    functionValue = cell2mat(cellfun(@(x) x.functionValue, population,'UniformOutput',false));
    constrViolation = cell2mat(cellfun(@(x) x.constrViolation, population,'UniformOutput',false));
    equalityConstrVals = cell2mat(cellfun(@(x) x.equalityConstrVals, population,'UniformOutput',false));
    inequalityConstrVals = cell2mat(cellfun(@(x) x.inequalityConstrVals, population,'UniformOutput',false));
    
    if parentSelectionType == 1
        [indicesParents, fitnessParents, functionValueParents, constrViolationParents, equalityConstrValsParents, inequalityConstrValsParents, parents]=selectParents(ulPop, fitnessVals, functionValue, constrViolation, equalityConstrVals, inequalityConstrVals, numParents, tag);
    elseif parentSelectionType == 2
        [indicesParents, fitnessParents, functionValueParents, constrViolationParents, equalityConstrValsParents, inequalityConstrValsParents, parents]=selectParents2(ulPop, fitnessVals, functionValue, constrViolation, equalityConstrVals, inequalityConstrVals, numParents, tag);
    elseif parentSelectionType == 3
        [indicesParents, fitnessParents, functionValueParents, constrViolationParents, equalityConstrValsParents, inequalityConstrValsParents, parents]=tournamentSelectParents(ulPop, fitnessVals, functionValue, constrViolation, equalityConstrVals, inequalityConstrVals, numParents, tag);
    elseif parentSelectionType == 4
        [indicesParents, fitnessParents, functionValueParents, constrViolationParents, equalityConstrValsParents, inequalityConstrValsParents, parents]=tournamentSelectSortedParents(ulPop, fitnessVals, functionValue, constrViolation, equalityConstrVals, inequalityConstrVals, numParents, tag);
    end

function [indicesParents, fitnessParents, functionValueParents, constrViolationParents, equalityConstrValsParents, inequalityConstrValsParents, parents]=selectParents(ulPop, fitnessVals, functionValue, constrViolation, equalityConstrVals, inequalityConstrVals, numParents, tag)
 
    %This function selects the best parent and other parents randomly for
    %crossover.
    ulPopSize = size(ulPop, 1);
        
    if sum(tag.ulPop==1)==0
        [~, index]=max(fitnessVals);
    else
        I = find(tag.ulPop==1);
        [~, index]=max(fitnessVals(tag.ulPop==1));
        index = I(index);
    end
    
    r = randint(1,numParents-1,[1 ulPopSize]);

    parents = [ulPop(index,:); ulPop(r,:)];
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
    
function [indicesParents, fitnessParents, functionValueParents, constrViolationParents, equalityConstrValsParents, inequalityConstrValsParents,parents]=selectParents2(ulPop, fitnessVals, functionValue, constrViolation, equalityConstrVals, inequalityConstrVals, numParents, tag)
 
    %This function selects the best parent and other parents by tournament for crossover.
    ulPopSize = size(ulPop, 1);
    numOtherParents = numParents-1;
    
    if sum(tag.ulPop==1)==0
        [~, index]=max(fitnessVals);
    else
        I = find(tag.ulPop==1);
        [~, index]=max(fitnessVals(tag.ulPop==1));
        index = I(index);
    end
    
    if ~isempty(constrViolation)
        t = find(constrViolation==0);
        if size(t,1)>2*numOtherParents
            permut = randperm(size(t,1));
            s = t(permut(1:2*numOtherParents));
        else
            [~, I] = sort(fitnessVals, 'descend');
            s = I(1:2*numOtherParents);
        end
    else
        permut = randperm(ulPopSize);
        s = permut(1:2*numOtherParents);
    end
    
    for i=1:2:2*numOtherParents
        if (fitnessVals(s(i))>fitnessVals(s(i+1)) && tag.ulPop(s(i))==tag.ulPop(s(i+1))) || (tag.ulPop(s(i))>tag.ulPop(s(i+1)))
            r((i+1)/2) = s(i);
        else
            r((i+1)/2) = s(i+1);
        end
    end

    parents = [ulPop(index,:); ulPop(r,:)];
    fitnessParents = [fitnessVals(index); fitnessVals(r)];
    functionValueParents = [functionValue(index,:); functionValue(r,:)];
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
    
function [indicesParents, fitnessParents, functionValueParents, constrViolationParents, equalityConstrValsParents, inequalityConstrValsParents, parents]=tournamentSelectParents(pop, fitnessVals, functionValue, constrViolation, equalityConstrVals, inequalityConstrVals, numParents, tag)
    
    %This function performs a tournament selection.
    %Member with tag=1 are preferred over members with tag=0.
    %If tag values are same then fitness comparison is performed.
    ulPopSize = size(pop, 1);
    
    if ~isempty(constrViolation)
        t = find(constrViolation==0);
        if size(t,1)>2*numParents
            permut = randperm(size(t,1));
            s = t(permut(1:2*numParents));
        else
            [~, I] = sort(fitnessVals, 'descend');
            s = I(1:2*numParents);
        end
    else
        permut = randperm(ulPopSize);
        s = permut(1:2*numParents);
    end
    
    for i=1:2:2*numParents
        if (fitnessVals(s(i))>fitnessVals(s(i+1)) && tag.ulPop(s(i))==tag.ulPop(s(i+1))) || (tag.ulPop(s(i))>tag.ulPop(s(i+1)))
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
    
function [indicesParents, fitnessParents, functionValueParents, constrViolationParents, equalityConstrValsParents, inequalityConstrValsParents, parents]=tournamentSelectSortedParents(pop, fitnessVals, functionValue, constrViolation, equalityConstrVals, inequalityConstrVals, numParents, tag)
    
    %This function performs a tournament selection.
    %Member with tag=1 are preferred over members with tag=0.
    %If tag values are same then fitness comparison is performed.
    %Parents are sorted by fitness before they are returned.
    ulPopSize = size(pop, 1);
    
    if ~isempty(constrViolation)
        t = find(constrViolation==0);
        if size(t,1)>2*numParents
            permut = randperm(size(t,1));
            s = t(permut(1:2*numParents));
        else
            [~, I] = sort(fitnessVals, 'descend');
            s = I(1:2*numParents);
        end
    else
        permut = randperm(ulPopSize);
        s = permut(1:2*numParents);
    end
    
    for i=1:2:2*numParents
        if (fitnessVals(s(i))>fitnessVals(s(i+1)) && tag.ulPop(s(i))==tag.ulPop(s(i+1))) || (tag.ulPop(s(i))>tag.ulPop(s(i+1)))
            r((i+1)/2) = s(i);
        else
            r((i+1)/2) = s(i+1);
        end
    end
        
    [~, I] = sort(fitnessVals(r),'descend');
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
    
function [offsprings] = geneticOperation (parents,probCrossover,probMutation,numOffsprings,crossoverType,mutationType)
    
    if crossoverType == 1 
        %planar pcx crossover
        offsprings=planar_pcx_crossover(parents, probCrossover, numOffsprings);
    elseif crossoverType == 2 
        %normal pcx crossover
        offsprings=pcx_crossover(parents, probCrossover, numOffsprings);
    else
        %Adaptive planar pcx crossover
        offsprings=adaptive_planar_pcx_crossover(parents, probCrossover, numOffsprings);
    end
    
    if mutationType == 1 
        %polynomial mutation
        offsprings=mutation(offsprings, probMutation);
    elseif mutationType == 2
        %variance based polynomial mutation
        offsprings=mutationVarianceBased(offsprings, probMutation, ulMemberVariance);       
    end   
    
function offsprings = planar_pcx_crossover(parents, probCrossover, numOffsprings)
 
    %Planar PCX Crossover
    %Enter parents as a matrix of size [dim X numParents]
    %The code has been written such that the first parent is chosen as the
    %index parent. Any number of parents can be given as input. However,
    %the other two parents are chosen randomly. Only three parents can be
    %used to perform a crossover.
 
    sigma1 = .1;
    sigma2 = .1;
    ulDim = size(parents,2);
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
    
function offsprings = pcx_crossover(parents, probCrossover, numOffsprings)
 
    %PCX Crossover
    %Enter parents as a matrix of size [dim X numParents]
    %The code has been written such that the first parent is chosen as the
    %index parent. Any number of parents can be given as input.
 
    sigma1 = 0.1;
    sigma2 = 0.1;
    ulDim = size(parents,2);
    numParents = size(parents,1);
    
    g = mean(parents);
 
    d_g = g - parents(1,:);
    d_g_sq = abs(d_g).^2;
    mod_d_g = sum(d_g_sq,2).^(1/2);
    
    d = parents(2:numParents,:) - parents(ones(1,numParents-1),:);
    d_sq = abs(d).^2;
    mod_d = sum(d_sq,2).^(1/2);
    
    if mod_d_g == 0 || all(mod_d) == 0
        for i=1:numOffsprings
            offsprings(i,:) = mutation(parents(1,:),1);
        end
    else
        theta = dot(d,d_g(ones(1,numParents-1),:),2)./(mod_d.*mod_d_g);
        theta(theta>1) = 1;
        rho = mean(mod_d.*sqrt(1-theta.^2));
        for i=1:numOffsprings
            rand_vect = rho*randn(1,ulDim)*sigma1;
            mod_rand_vect = norm(rand_vect);
 
            perp_vect = rand_vect - d_g*dot(rand_vect,d_g)/(mod_d_g*mod_d_g);
            para_vect = randn*sigma2*d_g;
            offsprings(i,:) = parents(1,:) + perp_vect + para_vect;
        end
    end
 
function offsprings = adaptive_planar_pcx_crossover(parents, probCrossover, numOffsprings)

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
    ulDim=size(offsprings,2);
    mum=20;
    for i=1:numOffsprings
        r = rand(1,ulDim);
        index = r<0.5;
        delta(index) = (2*r(index)).^(1/(mum+1)) - 1;
        index = r>=0.5;
        delta(index) = 1 - (2*(1 - r(index))).^(1/(mum+1));
 
        % Generate the corresponding child element.
        r = rand(1,ulDim);
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
    
function [population,fitnessVals,llOptimaFitness, tag] = update(population,children,llOptimaFitness, llOptimaFitnessOffsprings,tag,tagChildren)
    
    chooseMembers = max(size(children,1),2);
    if sum(tag.ulPop==1) <= 0.5*length(tag.ulPop)
        iTag0 = find(tag.ulPop == 0);
        permut = randperm(length(iTag0));
        r = permut(1:chooseMembers);
    else
        permut = randperm(size(population,1));
        r = permut(1:chooseMembers);
    end
    
    pool = [population(r); children];
    fitnessVals = cellfun(@(x) x.fitness, population);
    fitnessPool = cellfun(@(x) x.fitness, pool);
    tagPool = [tag.ulPop(r) tagChildren];
    llOptimaFitnessPool = [llOptimaFitness(r); llOptimaFitnessOffsprings];
    
    % Ensures that if maxIndex in the current (parent) population is selected, it does not get lost.
    % It also ensures that if something (children) better than maxIndex with tag=1 comes up, it does not get lost
    tempFitnessPool = fitnessPool;
    
    if sum(tag.ulPop==1)~=0    
        I = find(tag.ulPop==1);
        [~, index]=max(fitnessVals(tag.ulPop==1));
        maxIndex = I(index);
        if sum(maxIndex==r) == 1
            I = find(tagPool==1);
            [~, index]= max(fitnessPool(tagPool==1));       
            maxPoolIndex = I(index);
            tempFitnessPool(maxPoolIndex) = Inf;
        elseif sum(tagPool)>0
            if sum(fitnessVals(maxIndex)<=tempFitnessPool(tagPool==1))>0
                I = find(tagPool==1);
                [~, index]= max(fitnessPool(tagPool==1));
                maxPoolIndex = I(index);
                tempFitnessPool(maxPoolIndex) = Inf;
            end
        end
    end
    
    [~,I] = sort(tempFitnessPool,'descend');
    population(r) = pool(I(1:chooseMembers));
    tag.ulPop(r) = tagPool(I(1:chooseMembers));
    llOptimaFitness(r) = llOptimaFitnessPool(I(1:chooseMembers));
    fitnessVals(r) = fitnessPool(I(1:chooseMembers));


 
function offsprings=checkLimits(offsprings, ulDimMin, ulDimMax)
 
    numOffsprings = size(offsprings,1);
    dimMinMatrix = ulDimMin(ones(1,numOffsprings),:);
    offsprings(offsprings<dimMinMatrix)=dimMinMatrix(offsprings<dimMinMatrix);
    dimMaxMatrix = ulDimMax(ones(1,numOffsprings),:);
    offsprings(offsprings>dimMaxMatrix)=dimMaxMatrix(offsprings>dimMaxMatrix);
    
function [fitnessVals, constrViolation, population] = assignFitness(population, violationWindowEquality, violationWindowInequality, tag)

    functionValue = cell2mat(cellfun(@(x) x.functionValue, population,'UniformOutput',false));
    equalityConstrVals = cell2mat(cellfun(@(x) x.equalityConstrVals, population,'UniformOutput',false));
    inequalityConstrVals = cell2mat(cellfun(@(x) x.inequalityConstrVals, population,'UniformOutput',false));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Specify the number of constraints
    numEqualityConstr = size(equalityConstrVals,2);
    numInequalityConstr = size(inequalityConstrVals,2);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if numEqualityConstr == 0 && numInequalityConstr == 0
        constrViolation = [];
    end
    if numEqualityConstr == 0 && numInequalityConstr ~= 0
        inequalityConstrVals = inequalityConstrVals - violationWindowInequality;
        inequalityConstrVals(inequalityConstrVals<0)=0;
        constrViolation = sum(inequalityConstrVals,2);
    end
    if numEqualityConstr ~= 0 && numInequalityConstr == 0
        equalityConstrVals(abs(equalityConstrVals)<=violationWindowEquality) = 0;
        equalityConstrVals(abs(equalityConstrVals)>violationWindowEquality) = equalityConstrVals(abs(equalityConstrVals)>violationWindowEquality)-violationWindowEquality;
        constrViolation = sum(abs(equalityConstrVals),2);
    end
    if numEqualityConstr ~= 0 && numInequalityConstr ~= 0
        inequalityConstrVals = inequalityConstrVals - violationWindowInequality;
        inequalityConstrVals(inequalityConstrVals<0)=0;
        
        equalityConstrVals(abs(equalityConstrVals)<=violationWindowEquality) = 0;
        equalityConstrVals(abs(equalityConstrVals)>violationWindowEquality) = equalityConstrVals(abs(equalityConstrVals)>violationWindowEquality)-violationWindowEquality;
        constrViolation = sum(inequalityConstrVals,2)+ sum(abs(equalityConstrVals),2);
    end
    
    %Fitness assignment for single objective
    if size(functionValue,2) == 1
        if (isempty(constrViolation))
            fitnessVals = functionValue;
            for i = 1:length(population)
                population{i}.fitness = fitnessVals(i);
            end
            return;
        end

        if (constrViolation>0)
            fitnessVals = -constrViolation;
            for i = 1:length(population)
                population{i}.fitness = fitnessVals(i);
                population{i}.constrViolation = constrViolation(i);
            end
            return;
        end

        %Otherwise
        fitnessVals(constrViolation>0,1) = min(functionValue(constrViolation==0)) - constrViolation(constrViolation>0);
        fitnessVals(constrViolation==0,1) = functionValue(constrViolation==0);
        for i = 1:length(population)
            population{i}.fitness = fitnessVals(i);
            population{i}.constrViolation = constrViolation(i);
        end
        return;
    end
    %Fitness assignment for multi-objective
    if size(functionValue,2)>1
        if (isempty(constrViolation))
            %Following 3 lines modified on 15112014. Now tag 0 and tag 1 members are ranked separately.
            %Tag 0 fitness values are always worse than tag 1 fitness values.
            paretoRank = zeros(size(functionValue,1),1);
            try
                paretoRank(tag==1) = nonDominatedSorting(functionValue(tag==1,:));
            catch
                disp('hi');
            end
                
            paretoRank(tag==0) = max(paretoRank(tag==1)) + nonDominatedSorting(functionValue(tag==0,:));
            crowdingDistance = calculateCrowdingDistance(functionValue, paretoRank); 
            crowdingDistance = exp(-crowdingDistance);
            fitnessVals = 1+max(paretoRank) - paretoRank - crowdingDistance;
            for i = 1:length(population)
                population{i}.fitness = fitnessVals(i);
            end
            return;
        end

        if (constrViolation>0)
            fitnessVals = -constrViolation;
            for i = 1:length(population)
                population{i}.fitness = fitnessVals(i);
                population{i}.constrViolation = constrViolation(i);
            end
            return;
        end

        %Otherwise
        fitnessVals(constrViolation>0,1) = - constrViolation(constrViolation>0);
        indices = find(constrViolation==0);

        paretoRank(indices) = nonDominatedSorting(functionValue(indices,:));
        crowdingDistance(indices) = calculateCrowdingDistance(functionValue(indices,:), paretoRank(indices));
        crowdingDistance(indices) = exp(-crowdingDistance(indices));
        fitnessVals(indices,1) = 1+max(paretoRank(indices)) - paretoRank(indices) - crowdingDistance(indices);
        for i = 1:length(population)
            population{i}.fitness = fitnessVals(i);
            population{i}.constrViolation = constrViolation(i);
        end
    end
    
function tagAcceptLS = ifAcceptLS(initialIndivLS, eliteLS, violationWindowEquality, violationWindowInequality)

    eliteFunctionValue = initialIndivLS.functionValue; 
    eliteEqualityConstrVals = initialIndivLS.equalityConstrVals;
    eliteInequalityConstrVals = initialIndivLS.inequalityConstrVals;
    
    eliteFunctionValueLS = eliteLS.functionValue;
    eliteEqualityConstrValsLS = eliteLS.equalityConstrVals;
    eliteInequalityConstrValsLS = eliteLS.inequalityConstrVals;
    
    eliteFeasibilityLS=1;
    if ~isempty(eliteEqualityConstrVals)
        eliteFeasibilityLS = 0;
        if eliteEqualityConstrValsLS<=violationWindowEquality
            eliteFeasibilityLS = 1;
        end
    end
    if ~isempty(eliteInequalityConstrVals)
        eliteFeasibilityLS = 0;
        if eliteInequalityConstrValsLS<=violationWindowInequality
            eliteFeasibilityLS = 1;
        end
    end

    eliteFeasibility=1;
    if ~isempty(eliteEqualityConstrVals)
        eliteFeasibility = 0;
        if eliteEqualityConstrVals<=violationWindowEquality
            eliteFeasibility = 1;
        end
    end
    if ~isempty(eliteInequalityConstrVals)
        eliteFeasibility = 0;
        if eliteInequalityConstrVals<=violationWindowInequality
            eliteFeasibility = 1;
        end
    end

    if eliteFeasibilityLS>eliteFeasibility
        tagAcceptLS = 1;
    else
        tagAcceptLS = 0;
    end
    if eliteFeasibilityLS==1 && eliteFeasibility==1
        if eliteFunctionValueLS>eliteFunctionValue
            tagAcceptLS = 1;
        end
    end
% new additional in BLEAQ2
function [archive] = archiveUpdate(archive,newData,archiveSize,duplicatesCheck)
    
    % No duplicates check by Defualt
    if nargin < 4
        duplicatesCheck = 0;
    end
    
    % modification - 2016/08/19 
    % archive now is a cell array of structure
    archive = [archive; newData];
    
    % check for NaN function Value in the code 
    functionValue = cell2mat(cellfun(@(x) x.functionValue, archive,'UniformOutput',false));
    indexNaN = isnan(functionValue);
    if sum(indexNaN)>0
        error('NaN in the code. Check for error.');
    end
    
    % check and remove duplicates
    if duplicatesCheck == 1
        upper = cell2mat(cellfun(@(x) x.upper, archive, 'UniformOutput',false));
        [~,uniqueIndex] = unique(upper,'rows');
        archive = archive(uniqueIndex);
    end
    
    % updating mechanism - pick latest
    if length(archive) > archiveSize
        archive(1) = [];
    end 
    
function [Flag] = improvementsFlag(ulPopSize,gen,ulMaxGens,StoppingCriteria,archive,minSizePsiLS,minSizePhiLS)
    
    Flag.localSearch = 0; Flag.recheck = 0; 
    Flag.update = 0; Flag.replace = 0;
    Flag.dontDoPhi = 0;
    
    terminationFlag = (StoppingCriteria == 1) || (gen == ulMaxGens);
    frequency = ulPopSize;
    frequencyOffsetLS = frequency - 1;
    frequencyOffsetRC = floor(frequency/2 - 1);
    
    if length(archive.tag1) >= minSizePsiLS
        if terminationFlag == 1
            Flag.localSearch = 1; % Run EXACT local search
        elseif mod(gen+frequencyOffsetLS,frequency) == 0
            Flag.localSearch = 2; % Run default local search
        elseif mod(gen+frequencyOffsetRC,frequency) == 0
            Flag.recheck = 1;
        end
    end
    
    if length([archive.tag1;archive.tag0]) < minSizePhiLS
        Flag.dontDoPhi = 1;
    end