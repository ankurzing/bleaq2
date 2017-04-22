function [ulPop] = initialization(ulDim, ulDimMin, ulDimMax, llDimMin, llDimMax, ulPopSize, testProblemName)

    significantDistance = 1e-6;
    data = createData(ulDimMin, ulDimMax, llDimMin, llDimMax, testProblemName, []);
    if data.totalConstraints>0
        display('Ensuring feasibility for the constraints')
        constrModel = createModelConstraints(data);
    end
    
    for i=1:ulPopSize
        if data.totalConstraints>0 && i<=ceil(ulPopSize/2)
            %Ensuring relaxed feasibility
            weights(1) = rand();
            weights(2) = 1-rand();
            
            ulPop(i,:) = ensureFeasibility(ulDimMin, ulDimMax, llDimMin, llDimMax, data, constrModel, weights, []);
            for j=1:i-1
                if computeDistance(ulPop(i,:), ulPop(j,:))<=significantDistance
                    ulPop(i,:) = mutation(ulPop(i,:),1);
                    ulPop(i,:) = checkLimitsReflection(ulPop(i,:),ulDimMin,ulDimMax);
                    break;
                end
            end
        else
            %Normal initialization
            ulPop(i,:) = ulDimMin + rand(1, ulDim).*(ulDimMax-ulDimMin);
        end
    end

function data = createData(ulDimMin, ulDimMax, llDimMin, llDimMax, testProblemName, member)

    %If the input parameter 'member' is empty then population is created in the entire
    %box constraint set. If the input parameter 'member' is not empty then
    %population is created around the 'member'

    dimMin = [ulDimMin llDimMin];
    dimMax = [ulDimMax llDimMax];
    ulDim = size(ulDimMin,2);
    llDim = size(llDimMin,2);
    dim = ulDim + llDim;
    popSize = (dim+1)*(dim+2)/2+5*dim;
    if ~isempty(member)
        for i=1:popSize-1
            memberVariance = (dimMax - dimMin)/5;
            pop(i,:) = member + (sqrt(memberVariance)).*randn(size(member));
        end
        pop(popSize,:) = member;
        pop=checkLimitsReflection(pop, dimMin, dimMax);
    else
        pop=initialize(popSize,dim,dimMin,dimMax);
        member = pop(popSize,:);
    end
    
    for i=1:popSize
        [ulFunctionVals(i,:), ulEquality, ulInequality]=ulTestProblem(pop(i,1:ulDim), pop(i,ulDim+1:end), testProblemName);
        if ~isempty(ulEquality)
            ulEqualityConstrVals(i,:) = ulEquality;
        else
            ulEqualityConstrVals = [];
        end
        if ~isempty(ulInequality)
            ulInequalityConstrVals(i,:) = ulInequality;
        else
            ulInequalityConstrVals = [];
        end
        [llFunctionVals(i,:), llEquality, llInequality]=llTestProblem(pop(i,ulDim+1:end), testProblemName, pop(i,1:ulDim));
        if ~isempty(llEquality)
            llEqualityConstrVals(i,:) = llEquality;
        else
            llEqualityConstrVals = [];
        end
        if ~isempty(llInequality)
            llInequalityConstrVals(i,:) = llInequality;
        else
            llInequalityConstrVals = [];
        end
    end
    
    data.totalConstraints = size(ulEqualityConstrVals,2)+size(ulInequalityConstrVals,2)+size(llEqualityConstrVals,2)+size(llInequalityConstrVals,2);
    
    data.ul.pop = pop(:,1:ulDim);
    data.ll.pop = pop(:,ulDim+1:end);
    
    data.ul.functionVals = ulFunctionVals;
    data.ul.equalityConstrVals = ulEqualityConstrVals;
    data.ul.inequalityConstrVals = ulInequalityConstrVals;
    
    data.ll.functionVals = llFunctionVals;
    data.ll.equalityConstrVals = llEqualityConstrVals;
    data.ll.inequalityConstrVals = llInequalityConstrVals;
    
function constrModel = createModelConstraints(data)

    pop = [data.ul.pop data.ll.pop];
    
    equalityConstrVals = [data.ul.equalityConstrVals, data.ll.equalityConstrVals];
    inequalityConstrVals = [data.ul.inequalityConstrVals, data.ll.inequalityConstrVals];

    numEqualityConstr = size(equalityConstrVals,2);
    numInequalityConstr = size(inequalityConstrVals,2);
    
    if numEqualityConstr ~=0
        for i=1:numEqualityConstr
            approx.equalityConstr{i} = quadApprox(equalityConstrVals(:,i), pop);
        end
    else
        approx.equalityConstr = [];
    end

    if numInequalityConstr ~= 0
        for i=1:numInequalityConstr
            approx.inequalityConstr{i} = quadApprox(inequalityConstrVals(:,i), pop);
        end
    else
        approx.inequalityConstr = [];
    end
    
    constrModel = approx;
    
function [ulMember, llMember] = ensureFeasibility(ulDimMin, ulDimMax, llDimMin, llDimMax, data, constrModel, weights, member)

    %If input parameter 'member' is provided then the math program is
    %solved with the 'member' as starting point otherwise any random point
    %from data is chosen

    epsilonZero = 1e-6;
    
    dimMin = [ulDimMin llDimMin];
    dimMax = [ulDimMax llDimMax];
    ulDim = size(ulDimMin,2);
    llDim = size(llDimMin,2);
    dim = ulDim + llDim;
        
    %Copying constraint model to approx
    approx = constrModel;
    
    %If weights are not provided then function to be optimized is taken as
    %0, which means that only feasibility wrt constraints is needed
    %otherwise a new function with upper and lower level weights is
    %constructed
    if isempty(weights)
        functionVals = zeros(popSize,1);
    else
        functionVals = weights(1)*data.ul.functionVals + weights(2)*data.ll.functionVals;
    end
    
    pop = [data.ul.pop data.ll.pop];
    popSize = size(pop,1);
    if isempty(member)
        r = ceil(rand*popSize);
        member = pop(r,:);
    else
        functionVals = weights(1)*data.ul.functionVals + weights(2)*data.ll.functionVals;
    end
        
    %Adding function model to approx
    approx.function = quadApprox(functionVals, pop);

    options = optimset('Algorithm','interior-point','Display','off','TolX',1e-10,'TolFun',1e-10);
    
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
        [eliteIndiv] = linprog(-approx.function.linear,A_inequality,b_inequality',A_equality,b_equality',dimMin,dimMax,member,optionsLinprog);
        eliteIndiv = eliteIndiv';
    else
        [eliteIndiv] = fmincon(@(x) -approximatedFunction(x,approx.function),member,[],[],[],[],dimMin,dimMax,@(x) approximatedConstraints(x,approx.equalityConstr,approx.inequalityConstr),options);
    end
    ulMember = eliteIndiv(1:ulDim);
    llMember = eliteIndiv(ulDim+1:end);
    
    
function approxFunctionValue = approximatedFunction(pop, parameters)

    approxFunctionValue = parameters.constant + pop*parameters.linear + pop*parameters.sqmatrix*pop';

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
    
    
function pop=initialize(popSize,dim,dimMin,dimMax)

    for i=1:popSize
        pop(i,:) = dimMin + rand(1, dim).*(dimMax-dimMin);
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
    
    
function members = mutation(members, probMutation)
 
    numOffsprings=size(members,1);
    ulDim=size(members,2);
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
            members(i,r<probMutation) = members(i,r<probMutation) + delta(r<probMutation);
        end
    end
    
function members=checkLimitsReflection(members, dimMin, dimMax)
    %This function reflects an infeasible point into the variable bounds. If the
    %point lies far away, it assigns it a random position in the bounds.
    numOfMembers = size(members,1);
    dimMinMatrix = dimMin(ones(1,numOfMembers),:);
    dimMaxMatrix = dimMax(ones(1,numOfMembers),:);
    i = 0;
    while sum(sum(members<dimMinMatrix)) || sum(sum(members>dimMaxMatrix))
        I = members<dimMinMatrix-(dimMaxMatrix-dimMinMatrix);
        J = members>dimMaxMatrix+(dimMaxMatrix-dimMinMatrix);
        randI = rand(size(I));
        randJ = rand(size(J));
        members(I) = dimMinMatrix(I) + randI(I).*(dimMaxMatrix(I)-dimMinMatrix(I));
        members(J) = dimMinMatrix(J) + randJ(J).*(dimMaxMatrix(J)-dimMinMatrix(J));
        members(members<dimMinMatrix)=members(members<dimMinMatrix) + 2*(dimMinMatrix(members<dimMinMatrix)-members(members<dimMinMatrix));
        members(members>dimMaxMatrix)=members(members>dimMaxMatrix) + 2*(dimMaxMatrix(members>dimMaxMatrix)-members(members>dimMaxMatrix));
    end
    
